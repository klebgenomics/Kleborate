import os
import sys
import shutil
import subprocess
import pathlib
import tempfile
import urllib.request
import warnings
warnings.filterwarnings("ignore")


def install_dependencies():
    """Ensures the required libraries for BIGSdb_downloader are present.
    """
    required = {
        "requests": "requests",
        "requests_oauthlib": "requests-oauthlib",
        "rauth": "rauth",
        "bigsdb_downloader": "bigsdb-downloader",
    }
    missing_pip_names = []
    for module_name, pip_name in required.items():
        try:
            __import__(module_name)
        except ImportError:
            missing_pip_names.append(pip_name)

    if missing_pip_names:
        print(f"\n--- Installing missing dependencies ({', '.join(missing_pip_names)}) ---")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", *missing_pip_names])
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to pip install {missing_pip_names}: {e}")
            return False

    for module_name in required:
        try:
            __import__(module_name)
        except ImportError as e:
            print(f"[ERROR] {module_name} not importable after install: {e}")
            return False

    return True

_RAUTH_PATCH_SRC = '''
try:
    import rauth.session as _rauth_session

    _orig_parse_optional_params = _rauth_session.OAuth1Session._parse_optional_params

    def _patched_parse_optional_params(self, oauth_params, req_kwargs):
        if req_kwargs.get("params") is None:
            req_kwargs["params"] = {}
        if req_kwargs.get("data") is None:
            req_kwargs["data"] = {}
        return _orig_parse_optional_params(self, oauth_params, req_kwargs)

    _rauth_session.OAuth1Session._parse_optional_params = _patched_parse_optional_params
except ImportError:
    pass
'''


def make_rauth_patch_sitecustomize():
    patch_dir = pathlib.Path(tempfile.mkdtemp(prefix="rauth_patch_"))
    (patch_dir / "sitecustomize.py").write_text(_RAUTH_PATCH_SRC)
    return patch_dir


def run_bigsdb_downloader(args):
    runner_src = _RAUTH_PATCH_SRC + '''
import sys
from bigsdb_downloader.main import main
sys.exit(main())
'''
    cmd = [sys.executable, "-c", runner_src] + args
    return subprocess.run(cmd, check=True)


def get_paths():
    """locate Kleborate and MiST data directories."""
    try:
        import kleborate
        import mist

        k_path = pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__cgmlst' / 'data'
        m_path = pathlib.Path(mist.__file__).parent / 'resources' / 'pubmlst'

        return k_path, m_path
    except ImportError as e:
        print(f"Error: Missing dependency. {e}")
        sys.exit(1)


def patch_mist_resource(mist_res_dir):
    """download bigsdb downloader"""
    helper_file = mist_res_dir / "download_bigsdb.py"
    if not helper_file.exists():
        print(f"\n[FIX] MiST helper missing at {helper_file}. Downloading...")
        mist_res_dir.mkdir(parents=True, exist_ok=True)
        url = "https://raw.githubusercontent.com/B-02/MiST/master/mist/resources/pubmlst/download_bigsdb.py"
        try:
            urllib.request.urlretrieve(url, helper_file)
        except Exception as e:
            print(f"Failed to patch MiST: {e}")



def setup_bigsdb_credentials(token_path, key_name, site_name):
    """Checks for tokens and runs setup if they are missing."""
    exists = False
    if token_path.exists():
        if token_path.is_dir() and any(token_path.iterdir()):
            exists = True
        elif token_path.is_file():
            exists = True

    if exists:
        print(f"\n[INFO] Tokens found at {token_path}. Skipping authentication.")
        return

    if not install_dependencies():
        print("\n[ERROR] Could not install/import the dependencies required for bigsdb-downloader.")
        print("Please install it manually and re-run this script:")
        print("pip install bigsdb-downloader")
        sys.exit(1)

    print("\n--- Pasteur Authentication Setup ---")
    try:
        run_bigsdb_downloader([
            "--key_name", key_name,
            "--site", site_name,
            "--db", "pubmlst_klebsiella_seqdef",
            "--setup",
        ])
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] bigsdb-downloader setup failed: {e}")
        sys.exit(1)




def main():
    print("Kleborate cgMLST Database Setup")
    print("==========================================")

    if not shutil.which('mist'):
        print("Error: 'mist' not found.")
        sys.exit(1)

    # Use the current working directory
    cwd = pathlib.Path.cwd()

    # define the token dir
    token_base = cwd / ".bigsdb_tokens"
    token_check = token_base / "access_tokens"

    target_dir, mist_res_dir = get_paths()
    patch_mist_resource(mist_res_dir)

    if target_dir.exists():
        resp = input(f"Overwrite existing data at {target_dir}? (y/n): ").lower()
        if resp == 'y':
            shutil.rmtree(target_dir)
        else:
            print("Setup cancelled."); sys.exit(0)
    target_dir.mkdir(parents=True, exist_ok=True)

    print("\n1) Standard download\n2) Latest Pasteur (With authentication)")
    mode = input("Select (1 or 2): ")

    scheme_url = "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18"
    raw_download_path = target_dir / "kleb_scgmlst_s"
    index_path = target_dir / "kleb_scgmlst_s-index"

    download_env = None

    if mode == '2':
        setup_bigsdb_credentials(token_check, "Pasteur", "Pasteur")

        patch_dir = make_rauth_patch_sitecustomize()
        download_env = os.environ.copy()
        existing_pythonpath = download_env.get("PYTHONPATH", "")
        download_env["PYTHONPATH"] = (
            str(patch_dir) + (os.pathsep + existing_pythonpath if existing_pythonpath else "")
        )

        download_cmd = [
            "mist", "download", "--downloader", "bigsdb_auth",
            "--url", scheme_url, "--output", str(raw_download_path),
            "--include-profiles", "--dir-tokens", str(token_base),
            "--key-name", "Pasteur", "--site", "Pasteur"
        ]
    else:
        download_cmd = [
            "mist", "download", "--downloader", "bigsdb",
            "--url", scheme_url, "--output", str(raw_download_path), "--include-profiles"
        ]

    print(f"\n--- Downloading to {raw_download_path} ---")
    try:
        subprocess.run(download_cmd, check=True, env=download_env)
    except subprocess.CalledProcessError:
        print(f"\n[ERROR] Download failed. Check tokens in: {token_base}")
        sys.exit(1)

    print("\n--- Indexing ---")
    subprocess.run(["mist", "index", "--fasta-list", str(raw_download_path / "fasta_list.txt"),
                    "--profiles", str(raw_download_path / "profiles.tsv"),
                    "--output", str(index_path), "--threads", "8"], check=True)

    print(f"\nSuccess! Database ready at: {index_path}")

if __name__ == "__main__":
    main()

