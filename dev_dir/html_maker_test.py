import webbrowser
from pathlib import Path

mol_path = Path("/home/kato/10.git_repos/molgeom/dev_dir/h2o_360.xyz")
with open(mol_path, "r") as f:
    xyz_mol = f.read()

html_content = f"""<!DOCTYPE html>
<html>

<head>
    <title>3D Molecule Viewer</title>
    <script src="./js/jquery.min.js"></script>
    <script src="./js/3Dmol-min.js"></script>
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

tmp_html = Path("index.html")
with open(tmp_html, "w") as f:
    f.write(html_content)

# this is wsl
print(tmp_html.resolve().as_uri())  # file:///home/kato/10.git_repos/molgeom/index.html
# file:///home/kato/10.git_repos/molgeom/index.html   should be    file://wsl.localhost/Ubuntu/home/kato/10.git_repos/molgeom/index.html
webbrowser.get("wslview").open(tmp_html.resolve().as_uri())

# remove tmp_html
tmp_html.unlink()
