# mcpath/ui.py
import os, re, requests, yaml, subprocess
from IPython.display import display, clear_output
import ipywidgets as W

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

    # Only one visible at a time
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
        img_bytes = requests.get(url, timeout=20).content
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
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
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
    display(W.VBox(children, layout=W.Layout(width="auto")))

    # -------------------- handlers --------------------
    def on_clear(_):
        pdb_code.value = ""
        pdb_upload.value.clear()
        file_lbl.value = "No file chosen"
        chain_id.value = ""
        email.value = ""
        pred_type.value = "functional"
        # reset mode displays
        _sync_mode()
        # reset toggles & values
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
                # basic validation
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                # obtain PDB bytes & name
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # ---- Save PDB to SAVE_DIR ----
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # ---- Build mcpath_input.txt rows based on mode ----
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
                        email.value.strip() or "-"
                    ]

                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")

                # ---- Optional POST payload (if target_url set) ----
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
                    print(f"Submitting to {target_url} …")
                    r = requests.post(target_url, data=data, files=files, timeout=180)
                    print("HTTP", r.status_code)
                    try: print("JSON:", r.json())
                    except Exception: print("Response (≤800 chars):\n", r.text[:800])

                # ---- Run Octave after submit? ----
                if bool(cfg.get("run_octave_after_submit", False)):
                    SAVE_DIR_REAL = os.path.dirname(save_path)

                    # 1) Ensure pdbreader.m exists (parser shim or URL)
                    rp_parser = os.path.join(SAVE_DIR_REAL, "pdbreader.m")
                    pdbreader_url = (cfg.get("pdbreader_url") or "").strip()
                    try:
                        if pdbreader_url:
                            rb = requests.get(pdbreader_url, timeout=30).content
                            with open(rp_parser, "wb") as f: f.write(rb)
                            print(f"Fetched pdbreader.m from URL to: {rp_parser}")
                        else:
                            # use shim
                            with open(rp_parser, "w") as f: f.write(PDBREADER_SHIM_TEXT)
                            print(f"Using built-in pdbreader shim at: {rp_parser}")
                    except Exception as _e:
                        with open(rp_parser, "w") as f: f.write(PDBREADER_SHIM_TEXT)
                        print(f"(Fallback) Using built-in pdbreader shim at: {rp_parser}")

                    # 2) Ensure readpdb.m exists (MAIN analysis)
                    rp_readpdb = os.path.join(SAVE_DIR_REAL, "readpdb.m")
                    have_readpdb = False

                    # Try local repo copy first
                    try:
                        local_candidate = "/content/mcpath-colab/matlab/readpdb.m"
                        if os.path.exists(local_candidate):
                            with open(local_candidate, "rb") as src, open(rp_readpdb, "wb") as dst:
                                dst.write(src.read())
                            have_readpdb = True
                            print(f"Copied readpdb.m from repo to: {rp_readpdb}")
                    except Exception:
                        pass

                    # If not, try downloading from config
                    if not have_readpdb:
                        readpdb_url = (cfg.get("readpdb_url") or "").strip()
                        if readpdb_url:
                            try:
                                rb = requests.get(readpdb_url, timeout=30).content
                                with open(rp_readpdb, "wb") as f: f.write(rb)
                                have_readpdb = True
                                print(f"Fetched readpdb.m from URL to: {rp_readpdb}")
                            except Exception as e:
                                print("WARN: failed to fetch readpdb.m from readpdb_url:", e)

                    if not have_readpdb:
                        print("❌ readpdb.m not found. Please either:")
                        print("   - add it to your repo at mcpath-colab/matlab/readpdb.m, or")
                        print("   - set readpdb_url in config/defaults.yaml to a raw link of readpdb.m")
                        return  # stop here; nothing to run

                    # 3) Write runner and add paths
                    runner_path = os.path.join(SAVE_DIR_REAL, "run_mcpath.m")
                    # Optional extra paths from YAML
                    extra_paths = cfg.get("matlab_paths", []) or []
                    addpath_lines = ""
                    for p in extra_paths:
                        p_esc = p.replace("'", "''")
                        addpath_lines += f"addpath(genpath('{p_esc}'));\n"
                    # Always add /content and repo
                    addpath_lines += "addpath(genpath('/content')); addpath(genpath('/content/mcpath-colab'));\n"

                    runner = f"""
try
  % add paths
  {addpath_lines}
  % sanity check for readpdb
  if exist('readpdb','file') ~= 2
    error('readpdb.m not on path or missing');
  end

  fid = fopen('mcpath_input.txt','r'); assert(fid>0,'cannot open mcpath_input.txt');
  a = {{}}; line = fgetl(fid);
  while ischar(line), a{{end+1,1}} = line; line = fgetl(fid); end
  fclose(fid);
  % a{{2}} = pdb filename, a{{3}} = global chain ID
  readpdb(a{{2}}, a{{3}});
catch err
  fiderr = fopen('error','w');
  if fiderr>0
    fprintf(fiderr,'ERROR: %s\\n', err.message);
    fclose(fiderr);
  end
  exit(1);  % propagate failure to Python
end
exit(0);
"""
                    with open(runner_path, "w") as f:
                        f.write(runner)

                    # 4) Call Octave
                    cd_escaped = SAVE_DIR_REAL.replace("'", "''")
                    proc = subprocess.run(
                        ["octave", "-qf", "--eval", f"cd('{cd_escaped}'); run('run_mcpath.m');"],
                        text=True, capture_output=True
                    )

                    if proc.returncode != 0:
                        print("❌ Octave failed. Return code:", proc.returncode)
                        if proc.stdout:
                            print("Octave stdout (tail):\n", proc.stdout[-1000:])
                        if proc.stderr:
                            print("Octave stderr (tail):\n", proc.stderr[-1000:])
                        # If an error file exists, show its content
                        err_file = os.path.join(SAVE_DIR_REAL, "error")
                        if os.path.exists(err_file):
                            try:
                                print("\n--- error file ---")
                                with open(err_file, "r") as ef:
                                    print(ef.read()[:2000])
                            except Exception:
                                pass
                    else:
                        print("✅ Octave finished successfully.")


            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
