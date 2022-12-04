from geometry import Figure
from typing import AnyStr
from constants import BOUNDARY_LIMIT


def read_object_file(path: AnyStr):
    """
    Reads the object file provided
        input:
            path: path to the onject file which is given
        output:
            num_vertices: number of total vertives of the figure
            num_faces: number of faces of the figure
            vertices: idx to coordinate values of all vertices of the figure
            faces: idx of points making the faces in groups of 3
    """
    vertices = {}
    faces = []
    figure = None
    with open(path, encoding="utf8") as f:
        line_number = 0
        for line in f:
            content = line.split("\n")[0].split(",")
            # First line is special, it contains number of vertices and faces of the figure.
            if line_number == 0:
                num_vertices = int(content[0])
                num_faces = int(content[1])
                line_number = 1
            else:  # Get all vertices.
                if line_number <= num_vertices:
                    vertices[int(content[0])] = [
                        float(content[1]),
                        float(content[2]),
                        float(content[3]),
                    ]
                    line_number += 1
                else:  # Get all face values.
                    faces.append([int(content[0]), int(content[1]), int(content[2])])
        figure = Figure(num_vertices, num_faces)
        figure.construct_figure(vertices, faces)


def driver():
    read_object_file("object.txt")
