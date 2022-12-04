# Problem Statement
Add click and drag mouse functionality such that while the mouse button is pressed, movement of the mouse rotates the object thusly:
1. Horizontal movement of the mouse rotates the 3D object about the window’s Y-axis.
2. Vertical movement of the mouse rotates the 3D object about the window’s X-axis.
3. Diagonal movement of the mouse is decomposed into vertical and horizontal components and rotates the 3D object accordingly as above.
4. The point of the object nearest to the observer follows the mouse’s direction

## Dependencies
1. matplotlib

## To run
1. python main.py

## Parameters to change
Change the object.txt file to change the diagram.

## Flow of code
1. Read Object file
2. Plot initial diagram
3. repeat until ended (
    1. detect mouse click & drag
    2. find displacement vector
    3. calculate the change in point coordinates
    4. iterate over all the point that are connected and draw line
    5. display new image \
)
