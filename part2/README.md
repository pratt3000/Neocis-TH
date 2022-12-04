# Problem Statement
Create a separate program that contains all the functionality of part 1.
 Additionally, make each of the visible faces of the object a solid, opaque blue color. Make the color smoothly vary between #00005F (when the surface is viewed on edge, i.e. the normal of the surface makes a 90 degree angle to the Z-axis) and #0000FF (when the surface is viewed flat, i.e. orthogonal to the Z-axis) based on the angle with the Z-axis, such that the face is displayed similarly to how a shader would display it.


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
    4. iterate over all the faces and determine if the face has a normal vector facing outwards and towards the user. 
    5. Calculate slope of each normal, calculate color of surface and draw all such faces
    6. Display image
)
