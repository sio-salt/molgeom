import subprocess
import platform
import tempfile
import webbrowser
import time
import shutil
from pathlib import Path
from string import Template
from importlib.resources import files


def is_jupyter_notebook():
    """
    Check if the current Python environment is running in a Jupyter Notebook.

    Returns:
        bool: True if running in Jupyter Notebook, otherwise False.
    """
    try:
        from IPython import get_ipython

        if "IPKernelApp" in get_ipython().config:
            return True
    except (ImportError, AttributeError):
        pass
    return False


def is_wsl():
    return "microsoft" in platform.uname().release.lower()


def resolve_path(file_path: Path) -> str:
    """
    If WSL, converts Linux paths to Windows-style paths.
    """

    if is_wsl():
        try:
            # use wslpath command to convert Linux path to Windows path
            result = subprocess.run(
                ["wslpath", "-w", str(file_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
            win_path = result.stdout.strip()
            if win_path.startswith("\\\\"):
                win_path = win_path[2:]
            win_path = win_path.replace("\\", "/")
            return f"file://{win_path}"
        except subprocess.CalledProcessError as e:
            print(f"Error running wslpath: {e}")
            return file_path.as_uri()
    else:
        return file_path.resolve().as_uri()


def gen_mol_view_html(xyz_mol_data: str) -> str:
    """
    Generates HTML code for viewing molecular geometries using 3Dmol.js.
    args:
        xyz_mol_data: str
            XYZ-format molecular geometry data. can contain multiple molecules.
    returns:
        str: HTML string
    """

    html_template_path = files("molgeom").joinpath("templates/mol_view_template.html")
    with open(html_template_path, "r", encoding="utf-8") as f:
        html_template = Template(f.read())

    mol_js_uri = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"
    unique_id = str(time.time()).replace(".", "")
    return html_template.safe_substitute(
        {
            "mol_js_uri": mol_js_uri,
            "xyz_mol_data": xyz_mol_data,
            "unique_id": unique_id,
        }
    )


def open_html_in_browser(html_str: str, cleanup: bool = True) -> None:
    """
    Opens HTML string in a browser.
    args:
        html_str: str
            HTML string
        cleanup: bool
            If True, removes temporary file after opening (default: True)
    """

    tmp_dir = tempfile.mkdtemp()
    tmp_dir_path = Path(tmp_dir)
    # unique_id = str(time.time()).replace(".", "")
    # tmp_html_path = tmp_dir_path.joinpath(f"tmp_mol_view_{unique_id}.html")
    tmp_html_path = tmp_dir_path.joinpath("tmp_mol_view.html")
    tmp_html_path.write_text(html_str)

    if is_wsl():
        win_html_path = resolve_path(tmp_html_path)
        subprocess.run(["explorer.exe", win_html_path])
    else:
        uri = tmp_html_path.as_uri()
        webbrowser.open(uri)

    if cleanup:
        time.sleep(10)
        shutil.rmtree(tmp_dir)


def view_mol(
    xyz_mol_data: str, cleanup: bool = True, prefer_notebook: bool = True
) -> None:
    """
    View molecular geometry using 3Dmol.js in a browser.
    args:
        xyz_mol_data: str
            XYZ-format molecular geometry data. can contain multiple molecules.
        cleanup: bool
            If True, removes temporary file after opening (default: True)
        prefer_notebook: bool
            If True, tries to display in Jupyter Notebook if available (default: True)
    """

    display_in_notebook = prefer_notebook and is_jupyter_notebook()
    html_str = gen_mol_view_html(xyz_mol_data)

    if display_in_notebook:
        try:
            from IPython.display import display, HTML

            display(HTML(html_str))
        except (ImportError, ModuleNotFoundError):
            print("Failed to import and display in Jupyter Notebook.")
    else:
        open_html_in_browser(html_str, cleanup)


def test_view_mol():
    xyz_mol_data = """3
    H2O
    O  0.000000  0.000000  0.117370
    H  0.000000  0.755450 -0.469481
    H  0.000000 -0.755450 -0.469481
    5
    CH4
    C  0.000000  0.000000  0.000000
    H  0.000000 -0.000000  1.089000
    H  1.026719  0.000000 -0.363000
    H -0.513360 -0.889165 -0.363000
    H -0.513360  0.889165 -0.363000
    """
    view_mol(xyz_mol_data)


if __name__ == "__main__":
    test_view_mol()
