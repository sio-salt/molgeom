import os
import webbrowser

def generate_3dmol_html(molecule, output_path: str) -> None:
    """
    Generates an HTML file with 3Dmol.js to visualize the molecule.

    Args:
        molecule: Molecule instance containing atoms.
        output_path (str): Path to save the generated HTML file.
    """
    # Convert molecule data to JSON
    atom_data = [
        {"elem": atom.symbol, "x": atom.x, "y": atom.y, "z": atom.z}
        for atom in molecule
    ]
    
    # HTML template with embedded 3dmol.js
    threeDmol_js_url = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>3D Molecular Visualization</title>
        <script src="{threeDmol_js_url}"></script>
    </head>
    <body>
        <h1>Molecular Structure</h1>
        <div id="molViewer" style="width: 800px; height: 600px; position: relative;"></div>
        <script>
            // Initialize 3Dmol viewer
            const viewer = $3Dmol.createViewer("molViewer", {{ backgroundColor: "white" }});

            // Add atoms to the viewer
            const atoms = {atom_data};
            viewer.addAtoms(atoms);

            // Render the molecule
            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """

    # Write the HTML content to a file
    with open(output_path, "w") as file:
        file.write(html_content)
    print(f"HTML file generated at: {output_path}")

def open_html_file(filepath: str) -> None:
    """
    Open an HTML file in the default web browser.

    Args:
        filepath (str): Path to the HTML file.
    """
    webbrowser.open(f"file://{os.path.abspath(filepath)}")

# 使用例
if __name__ == "__main__":
    import subprocess
    from pathlib import Path
    from molgeom import Molecule
    
    # Create a sample molecule (replace with actual molecule data)
    
    # Generate the HTML file
    git_root = Path(subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip())
    dev_dir = git_root / "dev_dir"
    output_file = dev_dir / "molecule_visualization.html"
    molecule = Molecule.from_file(git_root / "tests/files" / "H2O.xyz")
    generate_3dmol_html(molecule, output_file)

    # Open the file in the default browser
    open_html_file(output_file)

