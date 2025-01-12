import subprocess
from pathlib import Path
import platform
import tempfile
import webbrowser
import time
from string import Template
from importlib.resources import files


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
    html_template_path = files("molgeom").joinpath("template/mol_view_template.html")
    with open(html_template_path, "r") as f:
        html_template = Template(f.read())

    jquery_js_path = resolve_path(files("molgeom").joinpath("static/js/jquery-3.7.1.min.js"))
    mol_js_path = resolve_path(files("molgeom").joinpath("static/js/3Dmol-2.4.2.min.js"))
    return html_template.substitute(
        {"jquery_js_path": jquery_js_path, "mol_js_path": mol_js_path, "xyz_mol_data": xyz_mol_data}
    )


def open_html_in_browser(html_str: str) -> None:
    """
    Opens HTML string in a browser.
    args:
        html_str: str
            HTML string
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_html_path = Path(tmp_dir).joinpath("tmp_mol_view.html")
        tmp_html_path.write_text(html_str)
        if is_wsl():
            win_html_path = resolve_path(tmp_html_path)
            subprocess.run(["explorer.exe", win_html_path])
        else:
            uri = tmp_html_path.as_uri()
            webbrowser.open(uri)

        time.sleep(3)


def view_mol(xyz_mol_data: str) -> None:
    """
    View molecular geometry using 3Dmol.js in a browser.
    args:
        xyz_mol_data: str
            XYZ-format molecular geometry data. can contain multiple molecules.
    """
    html_str = gen_mol_view_html(xyz_mol_data)
    open_html_in_browser(html_str)
