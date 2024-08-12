Users can load 3D models, rotate and scale models, and select two points to calculate the geodesic distance through mouse interaction on the GUI. 

The system development environment is VS2022.

Getting Started
Loading a 3D Model:
When you start the program, a window will open displaying a list of .obj files in the models folder.
Select a file from the list to load it into the program. The model will be displayed in the main view.

Interacting with the Model
 Selecting Points for Geodesic Distance Calculation
Left Mouse Button Click:
Click on the model to select points.
First Click: Select the start point (startPoint) for the geodesic distance calculation.
Second Click: Select the end point (endPoint).
After selecting both points, the program will automatically calculate and display the geodesic distance between them.

Rotating the Model
Right Mouse Button Drag:
Hold down the right mouse button and move the mouse to rotate the model.
Horizontal Movement: Rotates the model around the Y-axis (left-right rotation).
Vertical Movement: Rotates the model around the X-axis (up-down rotation).

Zooming and Panning the Camera
Zooming:
W Key: Moves the camera forward, zooming in on the model.
S Key: Moves the camera backward, zooming out from the model.
Panning:
A Key: Pans the camera to the left.
D Key: Pans the camera to the right.

Viewing the Results
The selected start and end points will be displayed in the ImGui interface along with their coordinates.
If both points are selected, the calculated geodesic distance will be displayed below the point coordinates.
