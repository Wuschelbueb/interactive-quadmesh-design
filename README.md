# Interactive Quad Mesh Design

This project is an advanced plugin designed specifically to integrate seamlessly with [OpenFlipper](https://gitlab.vci.rwth-aachen.de:9000/OpenFlipper-Free/OpenFlipper-Free.git), a powerful 3D modeling and animation tool. To utilize its capabilities, it is necessary to also include the [CoMISo](https://gitlab.vci.rwth-aachen.de:9000/CoMISo/CoMISo.git) library, which provides advanced computational geometry functionality. Together, this plugin and the CoMISo library provide a comprehensive set of tools to create Quad Meshes.

## Getting started

1. Create a new plugin folder in the OpenFlipper project
2. Clone the repository into the new folder
```bash
git clone https://github.com/Wuschelbueb/interactive-quadmesh-design.git
```
3. Make sure CoMISo is part of OpenFlipper
4. Build OpenFlipper with an IDE or editor

## How to use

It is extremely easy to use. First make sure that you loaded a 3D object of your choice.

1. Click the select button
    - Click once on a vertex to select it
    - Double-click anywhere to set the Cross Field direction
    - Set the distance - size of selection
    - Set the h parameter - size of the quads
2. Click the preview button for a preview
   - if you are not happy with the preview go back to step 1
3. Let the CoMISo solver solve the Mixed-Integer equations.
   - You end up with a Quad Mesh preview of the selection.
   - If you aren't happy with it, start at step 1 again.
   - If you would like to save it, go to step 4.
4. Save the 3D object with the texture anywhere you like.
5. Load it up with the [OpenGL](https://github.com/Wuschelbueb/opengl_visualization.git) visualization tool.

## License
>The MIT License (MIT)

>Copyright (c) 2023 wuschelbueb

>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
