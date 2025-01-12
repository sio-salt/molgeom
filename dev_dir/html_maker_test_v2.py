import subprocess
from pathlib import Path
import platform
import tempfile
import webbrowser
import shutil
import time


def create_temp_html(content, temp_dir):
    tmp_file = temp_dir / "temp_index.html"
    tmp_file.write_text(content)
    return tmp_file


def resolve_windows_path(file_path):
    """
    WSL環境下でLinuxパスをWindows形式のパスに変換し、file:// URIを生成する。
    非WSL環境では通常のfile:// URIを返す。
    """
    if "microsoft-standard" in platform.uname().release:
        try:
            # wslpathを使用してWindows形式のパスを取得
            result = subprocess.run(
                ["wslpath", "-w", str(file_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
            win_path = result.stdout.strip()  # 例: '\\wsl.localhost\Ubuntu\tmp\tmpabc123.html'
            if win_path.startswith("\\\\"):
                win_path = win_path[2:]  # 先頭の '\\' を削除
            win_path = win_path.replace("\\", "/")  # バックスラッシュをスラッシュに変換
            return f"file:///{win_path}"
        except subprocess.CalledProcessError as e:
            print(f"Error running wslpath: {e}")
            return file_path.as_uri()
    else:
        return file_path.resolve().as_uri()


def main():
    # スクリプトのディレクトリを取得
    script_dir = Path(__file__).parent.resolve()

    # jsディレクトリのパスを変数として定義
    js_dir = script_dir / "js"

    # jsディレクトリが存在するか確認
    if not js_dir.exists() or not js_dir.is_dir():
        print(f"JavaScript directory not found: {js_dir}")
        return

    # 必要なjsファイルのパスを定義
    jquery_js_filename = "jquery-3.7.1.min.js"
    jquery_js = js_dir / jquery_js_filename
    mol_js_filename = "3Dmol-2.4.2.min.js"
    mol_js = js_dir / mol_js_filename

    # jsファイルが存在するか確認
    if not jquery_js.exists() or not mol_js.exists():
        print(f"Required JavaScript files not found in {js_dir}")
        return

    # サンプルXYZファイルのパス
    mol_path = Path("/home/kato/10.git_repos/molgeom/dev_dir/h2o_360.xyz")

    # XYZファイルの内容を読み込む
    if mol_path.exists():
        with open(mol_path, "r") as f:
            xyz_mol = f.read()
    else:
        xyz_mol = "<!-- File not found -->"

    # 一時ディレクトリを作成
    temp_dir = Path(tempfile.mkdtemp())

    try:
        # jsディレクトリを一時ディレクトリにコピー
        temp_js_dir = temp_dir / "js"
        shutil.copytree(js_dir, temp_js_dir)

        # HTML内で使用するjsパスを変数として定義（相対パス）
        js_dir_relative = "./js"

        # HTMLコンテンツを定義
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>3D Molecule Viewer</title>
    <script src="{js_dir_relative}/{jquery_js_filename}"></script>
    <script src="{js_dir_relative}/{mol_js_filename}"></script>
    <style>
        .container {{
            display: flex;
            flex-direction: column;
            max-width: 800px;
            margin: 0 auto;
            padding: 20px;
            font-family: Arial, sans-serif;
        }}

        .viewer {{
            width: 100%;
            height: 500px;
            position: relative;
            box-sizing: border-box;
            border: 1px solid rgba(0, 0, 0, 0.15);
            border-radius: 8px;
            margin-bottom: 20px;
            background-clip: padding-box;
            overflow: hidden;
        }}

        .controls {{
            margin-bottom: 20px;
            display: flex;
            flex-direction: column;
            gap: 15px;
            padding: 15px;
            background: rgb(240, 240, 240);
            border-radius: 8px;
        }}

        .control-row {{
            display: flex;
            align-items: center;
            gap: 15px;
        }}

        select {{
            padding: 8px;
            border-radius: 4px;
            border: 1px solid #ccc;
            margin-right: 20px;
        }}

        input {{
            padding: 8px;
            border-radius: 4px;
            border: 1px solid #ccc;
        }}

        label {{
            font-weight: bold;
            color: black;
        }}

        #frameLabel {{
            font-family: monospace;
            min-width: 7ch;
            text-align: right;
            color: black;
        }}

        .nav-buttons {{
            display: flex;
            gap: 5px;
            margin-left: 10px;
        }}

        .nav-button,
        .toggle-button {{
            padding: 5px 10px;
            border: 1px solid #ccc;
            border-radius: 4px;
            background: white;
            cursor: pointer;
            font-weight: bold;
        }}

        .nav-button:hover,
        .toggle-button:hover {{
            background: #e0e0e0;
        }}

        .toggle-button.active {{
            background: #4CAF50;
            color: white;
        }}

        .toggle-button {{
            margin-left: 15px;
        }}
    </style>
</head>

<body>
    <div class="container">
        <div class="controls">
            <div class="control-row">
                <label for="frameSlider">Frame:</label>
                <input type="range" id="frameSlider" min="0" value="0" step="1" style="width: 200px; padding:0;">
                <span id="frameLabel"></span>
                <div class="nav-buttons">
                    <button class="nav-button" id="prevFrame">←</button>
                    <button class="nav-button" id="nextFrame">→</button>
                </div>
                <button id="toggleAxis" class="toggle-button">Show Axis</button>
            </div>
            <div class="control-row">
                <label for="styleSelect">Style:</label>
                <select id="styleSelect">
                    <option value="ball and stick">Ball and Stick</option>
                    <option value="stick">Stick</option>
                    <option value="sphere">Sphere</option>
                    <option value="line">Line</option>
                </select>
                <label for="bgColor">Background:</label>
                <select id="bgColor">
                    <option value="#333333">Dark Gray</option>
                    <option value="white">White</option>
                    <option value="black">Black</option>
                </select>
            </div>
        </div>
        <div id="viewer" class="viewer"></div>
    </div>

    <script>
        const molecule_data = `{xyz_mol}`;

        // Initialize viewer and model
        const viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
        const model = viewer.addModelsAsFrames(molecule_data, "xyz");

        // Get frame count and setup slider
        const frameCount = viewer.getNumFrames();
        const frameSlider = document.getElementById('frameSlider');
        frameSlider.max = frameCount - 1;

        let axisVisible = false;
        let currentLabels = [];

        function addAxis() {{
            const axisConfig = {{
                length: 5,
                radius: 0.075,
                bgColor: 'rgb(230, 230, 230)',
                bgOpacity: 0.95,
                colors: {{
                    x: 'rgb(255, 120, 120)',
                    y: 'rgb(120, 255, 120)',
                    z: 'rgb(120, 120, 255)'
                }},
                fontColors: {{
                    x: 'rgb(220, 60, 60)',
                    y: 'rgb(50, 200, 50)',
                    z: 'rgb(60, 60, 230)'
                }}
            }};

            const axes = ['x', 'y', 'z'];
            const positions = {{
                x: {{ start: {{ x: 0, y: 0, z: 0 }}, end: {{ x: axisConfig.length, y: 0, z: 0 }}, label: {{ x: axisConfig.length + 0.5, y: 0, z: 0 }} }},
                y: {{ start: {{ x: 0, y: 0, z: 0 }}, end: {{ x: 0, y: axisConfig.length, z: 0 }}, label: {{ x: 0, y: axisConfig.length + 0.5, z: 0 }} }},
                z: {{ start: {{ x: 0, y: 0, z: 0 }}, end: {{ x: 0, y: 0, z: axisConfig.length }}, label: {{ x: 0, y: 0, z: axisConfig.length + 0.5 }} }}
            }};

            axes.forEach(axis => {{
                viewer.addCylinder({{
                    start: positions[axis].start,
                    end: positions[axis].end,
                    radius: axisConfig.radius,
                    fromCap: true,
                    toCap: true,
                    color: axisConfig.colors[axis]
                }});
                viewer.addLabel(axis.toUpperCase(), {{
                    position: positions[axis].label,
                    backgroundColor: axisConfig.bgColor,
                    backgroundOpacity: axisConfig.bgOpacity,
                    fontColor: axisConfig.fontColors[axis],
                    borderThickness: 0.5,
                    borderColor: "black",
                    font: 'Arial',
                    fontSize: 20,
                }});
            }});

            // add origin
            viewer.addSphere({{
                center: {{x: 0 ,y: 0 ,z: 0 }},
                radius: 0.1,
                color: 'red'
            }});
        }}

        function removeAxis() {{
            viewer.removeAllShapes();
            viewer.removeAllLabels();
        }}

        function clearAtomLabels() {{
            currentLabels.forEach(label => viewer.removeLabel(label));
            currentLabels = [];
        }}

        // Event handlers for new controls
        document.getElementById('toggleAxis').addEventListener('click', function () {{
            axisVisible = !axisVisible;
            if (axisVisible) {{
                addAxis();
                this.classList.add('active');
            }} else {{
                removeAxis();
                this.classList.remove('active');
            }}
            viewer.render();
        }});

        // Previous event handlers and functions remain the same
        frameSlider.addEventListener('input', e => {{
            updateFrame(parseInt(e.target.value));
        }});

        document.getElementById('styleSelect').addEventListener('change', e => {{
            updateStyle(e.target.value);
        }});

        document.getElementById('bgColor').addEventListener('change', e => {{
            updateBgColor(e.target.value);
        }});

        function updateStyle(style) {{
            const styleConfig = style === "ball and stick"
                ? {{ stick: {{ radius: 0.125 }}, sphere: {{ scale: 0.22 }} }}
                : {{ [style]: {{}} }};

            model.setStyle({{}}, styleConfig);
            viewer.render();
        }}

        function updateBgColor(color) {{
            viewer.setBackgroundColor(color);
            viewer.render();
        }}

        const FRAME_SWITCH_INTERVAL = 80;
        const LONG_PRESS_DELAY = 500;
        let frameIntervalId = null;
        let longPressTimeoutId = null;

        function updateFrame(frame) {{
            frame = (frame + frameCount) % frameCount;
            frameSlider.value = frame;

            viewer.setFrame(frame)
            viewer.render();
            updateFrameLabel(frame);
        }}

        function startFrameChange(direction) {{
            if (frameIntervalId) return;

            const currentFrame = parseInt(frameSlider.value);
            updateFrame(currentFrame + direction);

            frameIntervalId = setInterval(() => {{
                const currentFrame = parseInt(frameSlider.value);
                updateFrame(currentFrame + direction);
            }}, FRAME_SWITCH_INTERVAL);
        }}

        function stopFrameChange() {{
            if (frameIntervalId) {{
                clearInterval(frameIntervalId);
                frameIntervalId = null;
            }}
            if (longPressTimeoutId) {{
                clearTimeout(longPressTimeoutId);
                longPressTimeoutId = null;
            }}
        }}

        function updateFrameLabel(frame) {{
            document.getElementById('frameLabel').textContent =
                `${{frame + 1}}/${{frameCount}}`;
        }}

        document.getElementById('prevFrame').addEventListener('mousedown', function () {{
            longPressTimeoutId = setTimeout(() => startFrameChange(-1), LONG_PRESS_DELAY);
            const currentFrame = parseInt(frameSlider.value);
            updateFrame(currentFrame - 1);
        }});

        document.getElementById('nextFrame').addEventListener('mousedown', function () {{
            longPressTimeoutId = setTimeout(() => startFrameChange(1), LONG_PRESS_DELAY);
            const currentFrame = parseInt(frameSlider.value);
            updateFrame(currentFrame + 1);
        }});

        ['mouseup', 'mouseleave'].forEach(eventName => {{
            document.getElementById('prevFrame').addEventListener(eventName, stopFrameChange);
            document.getElementById('nextFrame').addEventListener(eventName, stopFrameChange);
        }});

        ['touchstart', 'touchend', 'touchcancel'].forEach(eventName => {{
            document.getElementById('prevFrame').addEventListener(eventName, function (e) {{
                e.preventDefault();
                if (eventName === 'touchstart') {{
                    longPressTimeoutId = setTimeout(() => startFrameChange(-1), LONG_PRESS_DELAY);
                    const currentFrame = parseInt(frameSlider.value);
                    updateFrame(currentFrame - 1);
                }} else {{
                    stopFrameChange();
                }}
            }});

            document.getElementById('nextFrame').addEventListener(eventName, function (e) {{
                e.preventDefault();
                if (eventName === 'touchstart') {{
                    longPressTimeoutId = setTimeout(() => startFrameChange(1), LONG_PRESS_DELAY);
                    const currentFrame = parseInt(frameSlider.value);
                    updateFrame(currentFrame + 1);
                }} else {{
                    stopFrameChange();
                }}
            }});
        }});


        // Set initial states
        updateStyle("ball and stick");
        updateBgColor("#333333");
        viewer.zoomTo();
        viewer.render();
        updateFrameLabel(0);
    </script>
</body>

</html>
"""

        # 一時HTMLファイルを作成
        tmp_html = create_temp_html(html_content, temp_dir)

        # HTMLファイルをブラウザで開く
        if "microsoft-standard" in platform.uname().release:
            # WSL環境の場合
            win_html_path = resolve_windows_path(tmp_html).replace("file:///", "file://")
            print(f"Opening: {win_html_path}")
            try:
                subprocess.run(["explorer.exe", win_html_path])
            except FileNotFoundError:
                print("Error: explorer.exe not found.")
        else:
            # 通常の環境の場合
            uri = tmp_html.resolve().as_uri()
            print(f"Opening: {uri}")
            webbrowser.open(uri)

        # ブラウザがファイルを開くまで待機（必要に応じて調整）
        time.sleep(2)

    finally:
        try:
            shutil.rmtree(temp_dir)
            print(f"Deleted temporary directory: {temp_dir}")
        except Exception as e:
            print(f"Error deleting temporary directory: {e}")


if __name__ == "__main__":
    main()
