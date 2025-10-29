# mcpath/ui.py
import os, re, requests, yaml, subprocess, shutil, sys, glob, shlex
from IPython.display import display, clear_output
import ipywidgets as W

# ======================= SMART OCTAVE ENSURER (Colab + Drive cache + apt + micromamba) =======================
def _run(cmd, **kw):
    print("$", " ".join(cmd))
    return subprocess.run(cmd, check=True, text=True, **kw)

def _which_octave():
    return shutil.which("octave") or shutil.which("octave-cli")

def _mount_drive_if_colab():
    """Mount Google Drive if running on Colab and not mounted."""
    if not os.path.isdir("/content"):
        return False  # not Colab
    try:
        # If google.colab is importable, we are on Colab
        from google.colab import drive  # type: ignore
        # Check if already mounted
        if not os.path.isdir("/content/drive"):
            drive.mount("/content/drive", force_remount=False)
        else:
            # attempt a cheap access; if fails, re-mount
            try:
                _ = os.listdir("/content/drive")
            except Exception:
                drive.mount("/content/drive", force_remount=True)
        return True
    except Exception:
        return False

def _try_install_from_drive_cache(cache_dir="/content/drive/MyDrive/octave_debs"):
    """Install Octave from cached .deb files in Drive, if present."""
    if not os.path.isdir(cache_dir):
        print(f"Drive cache directory not found: {cache_dir}")
        return False
    debs = sorted(glob.glob(os.path.join(cache_dir, "*.deb")))
    if not debs:
        print(f"No .deb files in Drive cache: {cache_dir}")
        return False
    try:
        _run(["dpkg", "-i", *debs])
    except subprocess.CalledProcessError:
        # Resolve missing deps (may fetch a few packages)
        _run(["apt-get", "-f", "install", "-y", "-qq"])
    exe = _which_octave()
    if exe:
        print("✅ Octave installed from Drive cache:", exe)
        return True
    print("❌ Octave still not found after Drive cache install.")
    return False

def _build_or_refresh_drive_cache(cache_dir="/content/drive/MyDrive/octave_debs"):
    """Copy system APT .debs to Drive cache for future fast installs."""
    if not os.path.isdir("/content"):
        return  # non-Colab, skip
    os.makedirs(cache_dir, exist_ok=True)
    try:
        _run(["apt-get", "update", "-qq"])
        # Re-download only to ensure the cache has the exact versions
        _run(["apt-get", "install", "--reinstall", "-y", "-qq", "--download-only", "octave", "gnuplot"])
        # Copy .debs
        archives = "/var/cache/apt/archives"
        debs = [os.path.join(archives, f) for f in os.listdir(archives) if f.endswith(".deb")]
        if debs:
            _run(["bash", "-lc", f"cp -v {archives}/*.deb '{cache_dir}/'"])
            print(f"📦 Cached {len(debs)} .deb files to {cache_dir}")
    except Exception as e:
        print("WARN: could not refresh Drive cache:", e)

def _try_install_via_apt(install_timeout=300):
    """Install Octave using apt-get (Colab/Ubuntu)."""
    try:
        _run(["apt-get", "update", "-qq"])
        _run(["apt-get", "install", "-y", "-qq", "octave", "gnuplot"], timeout=install_timeout)
        exe = _which_octave()
        if exe:
            print("✅ Octave installed via apt:", exe)
            return True
    except Exception as e:
        print("❌ apt-get path failed:", e)
    return False

def _try_install_via_micromamba():
    """Micromamba fallback (creates env /usr/local/micromamba/envs/oct). Returns an executable shim if success."""
    try:
        print("⚙️ Setting up Octave via micromamba (conda-forge)…")
        mamba_root = "/usr/local/micromamba"
        env_name = "oct"
        mamba_bin = f"{mamba_root}/bin/micromamba"

        if not os.path.exists(mamba_bin):
            _run(["bash", "-lc",
                  "wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | "
                  "tar -xvj bin/micromamba >/dev/null 2>&1 || true"])
            os.makedirs(f"{mamba_root}/bin", exist_ok=True)
            _run(["bash", "-lc", f"mv -f bin/micromamba {mamba_root}/bin/"])

        env_prefix = f"{mamba_root}/envs/{env_name}"
        if not os.path.isdir(env_prefix):
            _run(["bash", "-lc", f"{mamba_bin} create -y -n {env_name} -c conda-forge octave gnuplot >/dev/null 2>&1 || true"])

        exe = shutil.which("octave") or shutil.which("octave-cli")
        # Even if system octave absent, we can run through micromamba
        shim = f"{mamba_bin} run -n {env_name} octave"
        print("✅ Micromamba octave shim ready.")
        return shim
    except Exception as e:
        print("❌ micromamba path failed:", e)
        return None

def ensure_octave_smart(
    drive_cache_dir="/content/drive/MyDrive/octave_debs",
    install_timeout=300,
    allow_micromamba=True
):
    """
    Smart installer:
      1) If octave exists, return it.
      2) If Colab: mount Drive.
      3) Try Drive cache dpkg install.
      4) If not installed, apt-get install; then refresh Drive cache.
      5) If still not installed (or not Colab), try micromamba.
      6) Return the path (or shim) to run octave. Raises if all fail.
    """
    exe = _which_octave()
    if exe:
        print("✅ Found Octave:", exe)
        return exe

    in_colab = os.path.isdir("/content")
    mounted = _mount_drive_if_colab() if in_colab else False

    # Try Drive cache
    if in_colab and mounted and _try_install_from_drive_cache(drive_cache_dir):
        exe = _which_octave()
        if exe:
            return exe

    # Try apt-get and then refresh cache for future speed
    if _try_install_via_apt(install_timeout=install_timeout):
        exe = _which_octave()
        if exe:
            if in_colab and mounted:
                _build_or_refresh_drive_cache(drive_cache_dir)
            return exe

    # Micromamba fallback
    if allow_micromamba:
        shim = _try_install_via_micromamba()
        if shim:
            return shim

    raise FileNotFoundError("GNU Octave not available. Tried PATH, Drive cache, apt-get, and micromamba.")

# -------------------- Octave pdbreader shim (fallback) --------------------
# Only used if you don't supply a real pdbreader.m via cfg['pdbreader_url'].
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
    def get_value(): return int(dd.value if mode.value == "list" else txt.value)
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
    show_title="MCPath-style Parameters"
):
    """Launch the MCPath parameter form."""
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

    # ---------- Mode 2 ----------
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
        minv=1, maxv=10_000_000, step=100, desc_style=DESC   # <-- maxv (not max)
    )

    # ---------- Mode 3 ----------
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                   description="Index of final residue:",
                                   layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:",
                         placeholder="B", layout=wide, style=DESC)
    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    (np_mode_3, np_container_3, np_dd_3, np_txt_3, get_num_paths_3) = _list_or_custom(
        label="Number of Paths", options=num_paths_opts_mode3, default_value=1000,
        minv=1, maxv=10_000_000, step=100, desc_style=DESC   # <-- maxv (not max)
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
    mode2_box = W.VBox([init_idx, init_chain, pl_mode_short, pl_container_short, np_mode_2, np_container_2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, np_mode_3, np_container_3])

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
    display(W.VBox(children, layout=W.Layout(width="auto")))

    # -------------------- handlers --------------------
    def on_clear(_):
        pdb_code.value = ""
        pdb_upload.value.clear()
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
                    rows = ["1", pdb_name, chain_global, str(get_big_len()), email.value.strip() or "-"]
                elif mode == "paths_init_len":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = ["2", pdb_name, chain_global, str(get_short_len()),
                            str(get_num_paths_2()), str(int(init_idx.value)),
                            (init_chain.value or "").strip(), email.value.strip() or "-"]
                else:
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = ["3", pdb_name, chain_global, str(int(init_idx.value)),
                            (init_chain.value or "").strip(), str(int(final_idx.value)),
                            (final_chain.value or "").strip(), email.value.strip() or "-"]

                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")

                # Optional POST
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip(): data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():    data[FN["email"]]    = email.value.strip()

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

                # ---- Run Octave after submit? (INSIDE on_submit) ----
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
                    except Exception as _e:
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

                    # Runner with diary + fail fast
                    runner_path = os.path.join(SAVE_DIR_REAL, "run_mcpath.m")
                    extra_paths = cfg.get("matlab_paths", []) or []
                    addpath_lines = ""
                    for p in extra_paths:
                        p_esc = p.replace("'", "''")
                        addpath_lines += f"addpath(genpath('{p_esc}'));\n"
                    addpath_lines += "addpath(genpath('/content')); addpath(genpath('/content/mcpath-colab'));\n"

                    runner = f"""
diary off; diary('octave.log');  % start logging
try
  {addpath_lines}
  if exist('readpdb','file') ~= 2, error('readpdb.m not on path'); end
  fid = fopen('mcpath_input.txt','r'); assert(fid>0,'cannot open mcpath_input.txt');
  a = {{}}; line = fgetl(fid);
  while ischar(line), a{{end+1,1}} = line; line = fgetl(fid); end
  fclose(fid);
  fprintf('Running readpdb on %%s (chain %%s)\\n', a{{2}}, a{{3}});
  readpdb(a{{2}}, a{{3}});
  diary off; exit(0);
catch err
  fiderr = fopen('error','w');
  if fiderr>0, fprintf(fiderr,'ERROR: %s\\n', err.message); fclose(fiderr); end
  diary off; exit(1);
end
"""
                    with open(runner_path, "w") as f:
                        f.write(runner)

                    # --- Ensure Octave (Drive cache → apt → micromamba) ---
                    print("▶ Launching Octave… (timeout=180s)")
                    try:
                        octave_cmd = ensure_octave_smart(
                            drive_cache_dir="/content/drive/MyDrive/octave_debs",
                            install_timeout=300,
                            allow_micromamba=True
                        )
                        cd_escaped = SAVE_DIR_REAL.replace("'", "''")
                        base_cmd = shlex.split(octave_cmd)  # supports micromamba shim
                        proc = subprocess.run(
                            base_cmd + ["-qf", "--eval", f"cd('{cd_escaped}'); run('run_mcpath.m');"],
                            text=True, timeout=180
                        )
                        rc = proc.returncode
                    except subprocess.TimeoutExpired:
                        print("❌ Octave timed out after 180s. Check /content/octave.log")
                        rc = -1
                    except FileNotFoundError as e:
                        print("❌", e)
                        rc = -2

                    # Show diary tail
                    log_path = os.path.join(SAVE_DIR_REAL, "octave.log")
                    if os.path.exists(log_path):
                        try:
                            with open(log_path, "r") as lf:
                                log = lf.read()
                            tail = log[-2000:] if len(log) > 2000 else log
                            print("\n--- octave.log (tail) ---\n", tail)
                        except Exception:
                            pass

                    if rc != 0:
                        err_file = os.path.join(SAVE_DIR_REAL, "error")
                        if os.path.exists(err_file):
                            with open(err_file, "r") as ef:
                                print("\n--- error file ---\n", ef.read()[:2000])
                        print("❌ Octave failed. Return code:", rc)
                    else:
                        print("✅ Octave finished successfully.")

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
