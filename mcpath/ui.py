# mcpath/ui.py
import os, re, requests, yaml, subprocess, shutil, sys
from IPython.display import display, clear_output
import ipywidgets as W

# -------------------- Singleton guard for the form --------------------
_FORM_SHOWN = False
_FORM_HANDLE = None

# -------------------- Octave availability helper --------------------
def _ensure_octave_and_pick_cmd(install_timeout=300):
    """
    Return a usable Octave executable.
    Resolution order:
      1) $OCTAVE_CMD env var
      2) PATH: octave, then octave-cli
      3) Common platform paths (Windows/macOS/Linux)
      4) Colab auto-install (apt-get) when running under /content
    Raises FileNotFoundError if nothing is found/installed.
    """
    import platform

    # 1) Allow explicit override
    override = os.environ.get("OCTAVE_CMD", "").strip()
    if override and shutil.which(override):
        return shutil.which(override)

    # 2) PATH lookup
    cand = shutil.which("octave") or shutil.which("octave-cli") or shutil.which("octave-cli.exe")
    if cand:
        return cand

    # 3) Common install locations by OS
    sysname = platform.system().lower()
    candidates = []

    if sysname == "windows":
        candidates += [
            r"C:\Octave\Octave-9.1.0\mingw64\bin\octave-cli.exe",
            r"C:\Octave\Octave-8.4.0\mingw64\bin\octave-cli.exe",
            r"C:\Program Files\GNU Octave\Octave-9.1.0\mingw64\bin\octave-cli.exe",
            r"C:\Program Files\GNU Octave\Octave-8.4.0\mingw64\bin\octave-cli.exe",
        ]
    elif sysname == "darwin":
        candidates += [
            "/opt/homebrew/bin/octave",        # Apple Silicon
            "/usr/local/bin/octave",           # Intel mac
            "/Applications/Octave.app/Contents/Resources/usr/bin/octave-cli",
        ]
    else:
        candidates += [
            "/usr/bin/octave", "/usr/local/bin/octave",
            "/usr/bin/octave-cli", "/usr/local/bin/octave-cli",
            os.path.expanduser("~/.conda/envs/base/bin/octave"),
            os.path.expanduser("~/.mambaforge/bin/octave"),
            os.path.expanduser("~/.miniforge3/bin/octave"),
        ]

    for p in candidates:
        if os.path.exists(p):
            return p

    # 4) Best-effort auto-install on Colab only
    in_colab = os.path.isdir("/content") and os.path.exists("/usr/bin/apt-get")
    if in_colab:
        try:
            print("⚙️ Octave not found — installing via apt-get (this may take a couple of minutes)…")
            subprocess.run(["apt-get", "update", "-qq"], check=True, text=True)
            subprocess.run(
                ["apt-get", "install", "-y", "-qq", "octave", "gnuplot"],
                check=True, text=True, timeout=install_timeout
            )
            cand = shutil.which("octave") or shutil.which("octave-cli")
            if cand:
                print("✅ Octave installed:", cand)
                return cand
        except Exception as e:
            print("❌ Failed to auto-install Octave on Colab:", e)

    # Nothing worked
    raise FileNotFoundError(
        "Octave executable not found.\n"
        "Install it and/or set OCTAVE_CMD to the full path.\n"
        "Quick install tips:\n"
        "  • Ubuntu/Debian:   sudo apt-get update && sudo apt-get install -y octave gnuplot\n"
        "  • macOS (Homebrew): brew install octave\n"
        "  • Windows (winget): winget install -e --id GNU.Octave   (or: choco install octave)\n"
        "  • Conda/Mamba:     conda install -c conda-forge octave gnuplot"
    )

# -------------------- Octave pdbreader shim (fallback) --------------------
PDBREADER_SHIM_TEXT = r"""
function pdb = pdbreader(filename)
  fid = fopen(filename,'r'); assert(fid>0,'pdbreader: cannot open %s',filename);
  Atoms = struct('AtomSerNo',{},'AtomName',{},'altLoc',{},'resName',{}, ...
                 'chainID',{},'resSeq',{},'X',{},'Y',{},'Z',{},'element',{},'iCode',{});
  k=0;
  while true
    t=fgetl(fid); if ~ischar(t), break; end
    if length(t)>=54 && (strncmp(t,'ATOM  ',6) || strncmp(t,'HETATM',6))
      k=k+1; s=@(a,b) strtrim(t(a:min(b,length(t))));
      Atoms(k).AtomSerNo=str2double(s(7,11));
      Atoms(k).AtomName =s(13,16);
      Atoms(k).altLoc   =s(17,17);
      Atoms(k).resName  =s(18,20);
      Atoms(k).chainID  =s(22,22);
      Atoms(k).resSeq   =str2double(s(23,26));
      Atoms(k).iCode    =s(27,27);
      Atoms(k).X        =str2double(s(31,38));
      Atoms(k).Y        =str2double(s(39,46));
      Atoms(k).Z        =str2double(s(47,54));
      if length(t)>=78
        elem=s(77,78);
      else
        nm=regexprep(Atoms(k).AtomName,'[^A-Za-z]',''); elem='';
        if ~isempty(nm), elem=nm(1); end
      end
      Atoms(k).element = elem;
    end
  end
  fclose(fid);
  pdb.Header = struct();
  pdb.Model  = struct('Atom', Atoms);
end
"""

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    print(f"Fetching PDB from RCSB: {url} (≤15s timeout)")
    r = requests.get(url, timeout=15); r.raise_for_status()
    return r.content

def _list_or_custom(label: str, options, default_value, minv, maxv, step=1, desc_style=None):
    """
    Switchable input: either a Dropdown (List) or a BoundedIntText (Custom).
    Returns (mode_toggle, container_box, dropdown_widget, int_text_widget, get_value_fn).
    """
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    common_layout = W.Layout(width="460px", min_width="420px")
    mode = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list",
        description="Input:",
        layout=W.Layout(width="260px"),
        style=desc_style
    )
    dd = W.Dropdown(
        options=options, value=default_value,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )
    txt = W.BoundedIntText(
        value=default_value, min=minv, max=maxv, step=step,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )

    container = W.VBox([dd])
    def _on_mode(change):
        container.children = [dd] if change["new"] == "list" else [txt]
    mode.observe(_on_mode, names="value")

    def get_value():
        return int(dd.value if mode.value == "list" else txt.value)

    return mode, container, dd, txt, get_value

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=10).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters",
    allow_multiple=False
):
    """Launch the MCPath parameter form."""
    global _FORM_SHOWN, _FORM_HANDLE

    # Singleton guard
    if _FORM_SHOWN and not allow_multiple:
        if _FORM_HANDLE is not None:
            display(_FORM_HANDLE)
        else:
            print("ℹ️ MCPath form already displayed.")
        return

    cfg = yaml.safe_load(requests.get(defaults_url, timeout=15).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    # Branding / labels
    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width': '260px'}
    wide = W.Layout(width="460px", min_width="420px")

    # ---------- PDB: code OR upload ----------
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")),
                        description="PDB code:",
                        layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose File")
    file_lbl   = W.Label("No file chosen")
    def _on_upload_change(change):
        if pdb_upload.value:
            fname = next(iter(pdb_upload.value.keys()))
            file_lbl.value = fname
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

    # ---------- Always-present fields ----------
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")),
                        description="Chain ID:",
                        layout=wide, style=DESC)
    email      = W.Text(value=str(cfg.get("email", "")),
                        description="Email (opt):",
                        layout=wide, style=DESC)

    # ---------- Prediction type ----------
    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues", "functional"),
            ("Allosteric Paths (initial residue + path length)", "paths_init_len"),
            ("Allosteric Paths (initial & final residues)", "paths_init_final"),
        ],
        value="functional",
        description="Prediction:",
        layout=W.Layout(width="auto"),
        style=DESC
    )

    # ---------- Functional mode: large path length ----------
    big_opts = cfg.get("path_length_options", [
        100000, 200000, 300000, 400000, 500000, 750000,
        1000000, 2000000, 3000000, 4000000
    ])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    (pl_mode_big, pl_container_big, pl_dd_big, pl_txt_big, get_big_len) = _list_or_custom(
        label="Path length", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000, desc_style=DESC
    )

    # ---------- Mode 2: initial residue + short path length + number of paths ----------
    init_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                  description="Index of initial residue:",
                                  layout=wide, style=DESC)
    init_chain = W.Text(value="", description="Chain of initial residue:",
                        placeholder="A", layout=wide, style=DESC)

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    (pl_mode_short, pl_container_short, pl_dd_short, pl_txt_short, get_short_len) = _list_or_custom(
        label="Length of Paths", options=short_len_opts, default_value=5,
        minv=1, maxv=10_000, step=1, desc_style=DESC
    )

    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000, 20000, 30000, 40000, 50000]
    (np_mode_2, np_container_2, np_dd_2, np_txt_2, get_num_paths_2) = _list_or_custom(
        label="Number of Paths", options=num_paths_opts_mode2, default_value=1000,
        minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    # ---------- Mode 3: initial & final residues + number of paths ----------
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                   description="Index of final residue:",
                                   layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:",
                         placeholder="B", layout=wide, style=DESC)

    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    (np_mode_3, np_container_3, np_dd_3, np_txt_3, get_num_paths_3) = _list_or_custom(
        label="Number of Paths", options=num_paths_opts_mode3, default_value=1000,
        minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    # ---------- Actions & output ----------
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    # ---------- Grouped layouts ----------
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start",
                                     flex_flow="row wrap", gap="10px"))

    functional_box = W.VBox([pl_mode_big, pl_container_big])
    mode2_box = W.VBox([
        init_idx, init_chain,
        pl_mode_short, pl_container_short,
        np_mode_2,   np_container_2
    ])
    mode3_box = W.VBox([
        init_idx, init_chain,
        final_idx, final_chain,
        np_mode_3,  np_container_3
    ])

    def _sync_mode(*_):
        functional_box.layout.display = "none"
        mode2_box.layout.display = "none"
        mode3_box.layout.display = "none"
        if pred_type.value == "functional":
            functional_box.layout.display = ""
        elif pred_type.value == "paths_init_len":
            mode2_box.layout.display = ""
        else:
            mode3_box.layout.display = ""
    _sync_mode()
    pred_type.observe(_sync_mode, names="value")

    # ---------- Render UI ----------
    children = []
    if logo: children.append(logo)
    children += [
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_row,
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box,
        mode2_box,
        mode3_box,
        chain_id, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]
    root = W.VBox(children, layout=W.Layout(width="auto"))
    display(root)
    _FORM_SHOWN = True
    _FORM_HANDLE = root

    # -------------------- handlers --------------------
    def on_clear(_):
        pdb_code.value = ""
        try:
            pdb_upload.value.clear()
        except Exception:
            pdb_upload.value = {}
        file_lbl.value = "No file chosen"
        chain_id.value = ""
        email.value = ""
        pred_type.value = "functional"
        _sync_mode()
        pl_mode_big.value = "list"
        pl_mode_short.value = "list"
        np_mode_2.value = "list"
        np_mode_3.value = "list"
        pl_dd_big.value = pl_dd_big.options[0]; pl_txt_big.value = int(pl_dd_big.options[0])
        init_idx.value = 1; init_chain.value = ""
        pl_dd_short.value = pl_dd_short.options[0]; pl_txt_short.value = int(pl_dd_short.options[0])
        np_dd_2.value = np_dd_2.options[0];       np_txt_2.value = int(np_dd_2.options[0])
        final_idx.value = 1; final_chain.value = ""
        np_dd_3.value = np_dd_3.options[0];       np_txt_3.value = int(np_dd_3.options[0])
        with out: clear_output()

    def _collect_pdb_bytes():
        if pdb_upload.value:
            (fname, meta) = next(iter(pdb_upload.value.items()))
            return meta["content"], fname
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_submit(_):
        with out:
            clear_output()
            try:
                print("▶ Validating inputs…")
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                print("▶ Getting PDB (upload or RCSB)…")
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # Save PDB
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # Build mcpath_input.txt
                input_path = os.path.join(os.path.dirname(save_path), "mcpath_input.txt")
                rows = []
                mode = pred_type.value

                if mode == "functional":
                    rows = [
                        "1",              # mode
                        pdb_name,         # pdb file name
                        chain_global,     # global chain
                        str(get_big_len()),
                        email.value.strip() or "-"
                    ]
                elif mode == "paths_init_len":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = [
                        "2",
                        pdb_name,
                        chain_global,
                        str(get_short_len()),
                        str(get_num_paths_2()),
                        str(int(init_idx.value)),
                        (init_chain.value or "").strip(),
                        email.value.strip() or "-"
                    ]
                else:  # paths_init_final
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = [
                        "3",
                        pdb_name,
                        chain_global,
                        str(int(init_idx.value)),
                        (init_chain.value or "").strip(),
                        str(int(final_idx.value)),
                        (final_chain.value or "").strip(),
                        email.value.strip() or "-",
                        str(int(get_num_paths_3())),
                    ]

                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")

                # Optional POST
                data = {
                    "prediction_mode": mode,
                    FN["chain_id"]: chain_global,
                }
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()

                if mode == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif mode == "paths_init_len":
                    data.update({
                        "length_paths":  int(get_short_len()),
                        "number_paths":  int(get_num_paths_2()),
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                    })
                else:
                    data.update({
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                        "index_final":   int(final_idx.value),
                        "chain_final":   (final_chain.value or "").strip(),
                        "number_paths":  int(get_num_paths_3()),
                    })

                files = {"pdb_file": (pdb_name, pdb_bytes, "chemical/x-pdb")}

                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data); preview["attached_file"] = pdb_name
                    print(preview)
                else:
                    print("▶ POSTing to server… (≤20s timeout)")
                    r = requests.post(target_url, data=data, files=files, timeout=20)
                    print("HTTP", r.status_code)
                    try: print("JSON:", r.json())
                    except Exception: print("Response (≤800 chars):\n", r.text[:800])

                # ---- Run Octave after submit? ----
                if bool(cfg.get("run_octave_after_submit", False)):
                    SAVE_DIR_REAL = os.path.dirname(save_path)

                    # Ensure pdbreader.m (shim or URL)
                    rp_parser = os.path.join(SAVE_DIR_REAL, "pdbreader.m")
                    pdbreader_url = (cfg.get("pdbreader_url") or "").strip()
                    try:
                        if pdbreader_url:
                            print(f"Fetching pdbreader.m (≤15s): {pdbreader_url}")
                            rb = requests.get(pdbreader_url, timeout=15).content
                            with open(rp_parser, "wb") as f: f.write(rb)
                            print(f"Saved: {rp_parser}")
                        else:
                            with open(rp_parser, "w") as f: f.write(PDBREADER_SHIM_TEXT)
                            print(f"Using built-in pdbreader shim at: {rp_parser}")
                    except Exception:
                        with open(rp_parser, "w") as f: f.write(PDBREADER_SHIM_TEXT)
                        print(f"(Fallback) Using built-in pdbreader shim at: {rp_parser}")

                    # Ensure readpdb.m
                    rp_readpdb = os.path.join(SAVE_DIR_REAL, "readpdb.m")
                    have_readpdb = False
                    try:
                        local_candidate = "/content/mcpath-colab/matlab/readpdb.m"
                        if os.path.exists(local_candidate):
                            with open(local_candidate, "rb") as src, open(rp_readpdb, "wb") as dst:
                                dst.write(src.read())
                            have_readpdb = True
                            print(f"Copied readpdb.m to: {rp_readpdb}")
                    except Exception:
                        pass
                    if not have_readpdb:
                        readpdb_url = (cfg.get("readpdb_url") or "").strip()
                        if readpdb_url:
                            try:
                                print(f"Fetching readpdb.m (≤15s): {readpdb_url}")
                                rb = requests.get(readpdb_url, timeout=15).content
                                with open(rp_readpdb, "wb") as f: f.write(rb)
                                have_readpdb = True
                                print(f"Saved: {rp_readpdb}")
                            except Exception as e:
                                print("WARN: failed to fetch readpdb.m:", e)
                    if not have_readpdb:
                        print("❌ readpdb.m not found. Set readpdb_url in config or place it in mcpath-colab/matlab/")
                        return

                    # === Ensure atomistic.m, infinite.m, and required data files ===
                    def _ensure_file(target_path, local_candidates=(), url_key=None, human_name=None):
                        """Write file to target_path from first available source: local_candidates then cfg[url_key]."""
                        os.makedirs(os.path.dirname(target_path), exist_ok=True)
                        # 1) local candidates
                        for lc in local_candidates:
                            try:
                                if lc and os.path.exists(lc):
                                    with open(lc, "rb") as src, open(target_path, "wb") as dst:
                                        dst.write(src.read())
                                    print(f"Saved {human_name or os.path.basename(target_path)} from local: {lc}")
                                    return True
                            except Exception as e:
                                print(f"WARN: copy failed for {lc}: {e}")
                        # 2) URL fallback
                        if url_key and cfg.get(url_key):
                            try:
                                url = cfg[url_key].strip()
                                print(f"Fetching {human_name or os.path.basename(target_path)} (≤15s): {url}")
                                rb = requests.get(url, timeout=15).content
                                with open(target_path, "wb") as f: f.write(rb)
                                print(f"Saved: {target_path}")
                                return True
                            except Exception as e:
                                print(f"WARN: failed to fetch {human_name or target_path}: {e}")
                        return False

                    atomistic_m   = os.path.join(SAVE_DIR_REAL, "atomistic.m")
                    infinite_m    = os.path.join(SAVE_DIR_REAL, "infinite.m")
                    vdw_param     = os.path.join(SAVE_DIR_REAL, "vdw_cns.param")
                    pdb_top       = os.path.join(SAVE_DIR_REAL, "pdb_cns.top")

                    _ensure_file(atomistic_m,
                                 local_candidates=["/content/mcpath-colab/matlab/atomistic.m"],
                                 url_key="atomistic_url", human_name="atomistic.m")
                    _ensure_file(infinite_m,
                                 local_candidates=["/content/mcpath-colab/matlab/infinite.m"],
                                 url_key="infinite_url", human_name="infinite.m")
                    _ensure_file(vdw_param,
                                 local_candidates=["/content/mcpath-colab/matlab/vdw_cns.param"],
                                 url_key="vdw_param_url", human_name="vdw_cns.param")
                    _ensure_file(pdb_top,
                                 local_candidates=["/content/mcpath-colab/matlab/pdb_cns.top"],
                                 url_key="pdb_top_url", human_name="pdb_cns.top")

                    for must in [atomistic_m, infinite_m, vdw_param, pdb_top]:
                        if not os.path.exists(must):
                            print(f"❌ Required file missing: {must}. Provide local copy or set its URL in config.")

                    # --------- Write runner (with visible progress) ----------
                    runner_path = os.path.join(SAVE_DIR_REAL, "run_mcpath.m")
                    extra_paths = cfg.get("matlab_paths", []) or []
                    addpath_lines = ""
                    for p in extra_paths:
                        p_esc = p.replace("'", "''")
                        addpath_lines += f"addpath(genpath('{p_esc}'));\n"
                    addpath_lines += "addpath(genpath('/content')); addpath(genpath('/content/mcpath-colab'));\n"

                    cd_escaped = SAVE_DIR_REAL.replace("'", "''")
                    runner = f"""
diary off; diary('octave.log');
try
  {addpath_lines}
  cd('{cd_escaped}');
  printf('RUNNER_START\\n'); fflush(stdout);

  if exist('readpdb','file') ~= 2, error('readpdb.m not on path'); end
  if exist('atomistic','file') ~= 2, warning('atomistic.m not found on path'); end
  if exist('infinite','file')  ~= 2, warning('infinite.m not found on path'); end

  fid = fopen('mcpath_input.txt','r'); assert(fid>0,'cannot open mcpath_input.txt');
  a = {{}}; line = fgetl(fid);
  while ischar(line), a{{end+1,1}} = strtrim(line); line = fgetl(fid); end
  fclose(fid);

  printf('RUNNER_GOT_INPUT: mode=%s pdb=%s chain=%s\\n', a{{1}}, a{{2}}, a{{3}}); fflush(stdout);

  fprintf('Running readpdb on %s (chain %s)\\n', a{{2}}, a{{3}});
  readpdb(a{{2}}, a{{3}});
  printf('READPDB_DONE\\n'); fflush(stdout);

  mode = str2double(a{{1}});
  steps = 0;
  if mode == 1
      steps = str2double(a{{4}});
  elseif mode == 2
      steps = str2double(a{{4}});
  elseif mode == 3
      if numel(a) >= 8, steps = str2double(a{{8}}); else, steps = 1000; end
  else
      steps = 1000;
  end
  if ~isfinite(steps) || steps <= 0, steps = 1000; end

  if exist('atomistic','file') == 2
      fprintf('Running atomistic(%s)\\n', a{{2}});
      try
          atomistic(a{{2}});
          printf('ATOMISTIC_DONE\\n'); fflush(stdout);
      catch atomErr
          warning('atomistic() failed: %s', atomErr.message);
      end
  else
      warning('atomistic.m missing — skipping atomistic()');
  end

  if exist('infinite','file') == 2
      fprintf('Running infinite(%s, %d)\\n', a{{2}}, steps);
      try
          infinite(a{{2}}, steps);
          printf('INFINITE_DONE\\n'); fflush(stdout);
      catch infErr
          warning('infinite() failed: %s', infErr.message);
      end
  else
      warning('infinite.m missing — skipping infinite()');
  end

  coor_src = [a{{2}} '.cor'];
  atom_src = [a{{2}} '_atomistic.out'];
  path_src = [a{{2}} '_atomistic_' num2str(steps) 'steps_infinite.path'];

  function tf = copy_to_free_name(src, dst_base)
      tf = false;
      if exist(src, 'file') ~= 2, return; end
      dst = dst_base; k = 1;
      while exist(dst, 'file') == 2
          dst = sprintf('%s.%d', dst_base, k); k = k + 1;
      end
      try, copyfile(src, dst); tf = true; fprintf('Saved %s -> %s\\n', src, dst);
      catch, warning('Failed to copy %s -> %s', src, dst);
      end
  end

  copy_to_free_name(coor_src, 'coor_file');
  copy_to_free_name(atom_src, 'atom_file');
  copy_to_free_name(path_src, 'path_file');

  printf('RUNNER_END\\n'); fflush(stdout);
  diary off; exit(0);
catch err
  fiderr = fopen('error','w');
  if fiderr>0, fprintf(fiderr,'ERROR: %s\\n', err.message); fclose(fiderr); end
  diary off; exit(1);
end
"""
                    with open(runner_path, "w") as f:
                        f.write(runner)

                    # --- Ensure Octave is present and run (hardened) ---
                    print("▶ Launching Octave… (timeout=420s)")
                    try:
                        octave_cmd = _ensure_octave_and_pick_cmd(install_timeout=420)
                        print(f"Using Octave at: {octave_cmd}")

                        # Version
                        try:
                            ver = subprocess.run([octave_cmd, "--version"], text=True, capture_output=True, timeout=30)
                            print("Octave version:\n", (ver.stdout or ver.stderr)[:2000])
                        except Exception as e_ver:
                            print("WARN: could not query Octave version:", e_ver)

                        # Smoke test
                        try:
                            smoke = subprocess.run(
                                [octave_cmd, "-qf", "--no-gui", "--no-window-system",
                                 "--eval", "disp('OCTAVE_SMOKE_OK'); fflush(stdout); exit(0);"],
                                text=True, capture_output=True, timeout=30
                            )
                            if smoke.returncode != 0 or "OCTAVE_SMOKE_OK" not in (smoke.stdout + smoke.stderr):
                                print("WARN: Octave smoke test output:\n", (smoke.stdout + smoke.stderr)[:1000])
                                raise RuntimeError("Octave smoke test failed")
                        except Exception as e_smoke:
                            print("❌ Octave smoke test failed:", e_smoke)
                            raise

                        # Run runner (correct cwd + self-contained eval)
                        eval_code = (
                            "try, "
                            f"  cd('{SAVE_DIR_REAL.replace(\"'\", \"''\")}'); "
                            "  disp('RUN_EVAL_CALLING_RUN_MCPATH'); fflush(stdout); "
                            "  run('run_mcpath.m'); "
                            "catch err, "
                            "  disp('RUNNER_ERROR_BEGIN'); "
                            "  disp(getReport(err, 'extended')); "
                            "  disp('RUNNER_ERROR_END'); "
                            "  exit(1); "
                            "end; "
                            "exit(0);"
                        )

                        proc = subprocess.run(
                            [octave_cmd, "-qf", "--no-gui", "--no-window-system", "--eval", eval_code],
                            text=True, capture_output=True, timeout=420,
                            cwd=SAVE_DIR_REAL
                        )
                        rc = proc.returncode
                        if proc.stdout:
                            print("\n--- Octave STDOUT (tail) ---\n", proc.stdout[-2000:])
                        if proc.stderr:
                            print("\n--- Octave STDERR (tail) ---\n", proc.stderr[-2000:])

                    except subprocess.TimeoutExpired:
                        print("❌ Octave timed out after 420s. Check /content/octave.log")
                        rc = -1
                    except FileNotFoundError as e:
                        print("❌", e)
                        rc = -2
                    except Exception as e:
                        print("❌ Unexpected error launching Octave:", e)
                        rc = -3

                    # Show diary tail and error file
                    log_path = os.path.join(SAVE_DIR_REAL, "octave.log")
                    if os.path.exists(log_path):
                        try:
                            with open(log_path, "r") as lf:
                                log = lf.read()
                            tail = log[-4000:] if len(log) > 4000 else log
                            print("\n--- octave.log (tail) ---\n", tail)
                        except Exception as e:
                            print("WARN: could not read octave.log:", e)

                    err_file = os.path.join(SAVE_DIR_REAL, "error")
                    if os.path.exists(err_file):
                        try:
                            with open(err_file, "r") as ef:
                                print("\n--- error file ---\n", ef.read()[:4000])
                        except Exception as e:
                            print("WARN: could not read error file:", e)

                    for pth in ["coor_file", "atom_file", "path_file"]:
                        print(f"Exists {pth}? ", os.path.exists(os.path.join(SAVE_DIR_REAL, pth)))

                    if rc != 0:
                        print("❌ Octave failed. Return code:", rc)
                    else:
                        print("✅ Octave finished successfully.")

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)

def close_form():
    """Close the current form so a new one can be created."""
    global _FORM_SHOWN, _FORM_HANDLE
    if _FORM_HANDLE is not None:
        try:
            _FORM_HANDLE.close()
        except Exception:
            pass
    _FORM_SHOWN = False
    _FORM_HANDLE = None

# --- Auto-launch in Colab so the form is ready to fill immediately ---
def _in_colab():
    try:
        import google.colab  # type: ignore
        return True
    except Exception:
        return False

def autolaunch(defaults_url: str = "https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml"):
    """Manually trigger the UI (useful outside Colab or if autolaunch is disabled)."""
    try:
        launch(defaults_url=defaults_url, allow_multiple=False)
    except Exception as e:
        print("❌ autolaunch() failed:", e)

# Auto-run when imported in Colab (can be disabled via MCPATH_AUTOLAUNCH=0)
if os.environ.get("MCPATH_AUTOLAUNCH", "1") == "1" and _in_colab():
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            print("🔧 Auto-launching MCPath form…")
            autolaunch()
    except Exception as _e:
        print("⚠️ Auto-launch skipped:", _e)
