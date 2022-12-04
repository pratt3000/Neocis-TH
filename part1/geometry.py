from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import math

from constants import (
    BOUNDARY_LIMIT,
    FIGURE_REFRESH_DELTA,
    TOTAL_ROTATION_PER_ITERATION,
)


class Point:
    def __init__(self, x_coordinate: float, y_coordinate: float, z_coordinate: float):
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate


class Face:
    def __init__(self, coordinates: Tuple[Point, Point, Point]):
        self.coordinates = coordinates


class Figure:
    def __init__(self, num_vertices, num_faces):
        self.num_vertices = num_vertices
        self.num_faces = num_faces
        self.vertices = {}
        self.faces = []
        self.press = None
        self.centroid_x = None
        self.centroid_y = None
        self.x_max = None
        self.y_max = None

    def _rotate_x_3D(self, theta: float):
        """
        Rotates all of the vertices of the given figure according to
        the theta given along the 'x-axis'
        Input:
            theta: The angle by which we need to rotate
        """
        sinTheta = math.sin(theta)
        cosTheta = math.cos(theta)

        # Rotating all of the vectors of depending on the theta given using basic projectio formula.
        for point_idx in self.vertices:
            point = self.vertices[point_idx]

            cur_x = point.x_coordinate
            cur_y = point.y_coordinate * cosTheta - point.z_coordinate * sinTheta
            cur_z = point.z_coordinate * cosTheta + point.y_coordinate * sinTheta

            self.vertices[point_idx] = Point(cur_x, cur_y, cur_z)

    def _rotate_y_3D(self, theta: float):
        """
        Rotates all of the vertices of the given figure according to the theta given along the 'y-axis'
        Input:
            theta: The angle by which we need to rotate
        """
        sinTheta = math.sin(theta)
        cosTheta = math.cos(theta)

        # Rotating all of the vectors of depending on the theta given using basic projectio formula.
        for point_idx in self.vertices:
            point = self.vertices[point_idx]

            cur_x = point.x_coordinate * cosTheta + point.z_coordinate * sinTheta
            cur_y = point.y_coordinate
            cur_z = point.z_coordinate * cosTheta - point.x_coordinate * sinTheta

            self.vertices[point_idx] = Point(cur_x, cur_y, cur_z)

    def _rotate_z_3D(self, theta: float):
        """
        ***
            We dont really need/use this function but this is there just for
            symmetry and incase its needed in the future
        ***
        Rotates all of the vertices of the given figure according to the theta given along the 'z-axis'
        Input:
            theta: The angle by which we need to rotate
        """
        sinTheta = math.sin(theta)
        cosTheta = math.cos(theta)

        # Rotating all of the vectors of depending on the theta given using basic projection formula.
        for point_idx in self.vertices:
            point = self.vertices[point_idx]

            cur_x = point.x_coordinate * cosTheta - point.y_coordinate * sinTheta
            cur_y = point.y_coordinate * cosTheta + point.x_coordinate * sinTheta
            cur_z = point.z_coordinate

            self.vertices[point_idx] = Point(cur_x, cur_y, cur_z)

    def _xy_max(self):
        """
        Calculate centroids and max values around a and y axis.
        """
        x_coordinates = [vertex.x_coordinate for vertex in self.vertices.values()]
        y_coordinates = [vertex.y_coordinate for vertex in self.vertices.values()]

        self.centroid_x = sum(x_coordinates) / len(x_coordinates)
        self.centroid_y = sum(y_coordinates) / len(y_coordinates)

        self.x_max = max([abs(x) for x in x_coordinates])
        self.y_max = max([abs(y) for y in y_coordinates])

    def _xy_boundary_limits(self):
        """
        Set the boundary limits of the bounding box of the displaying plot box.
        """
        plt.xlim(
            [
                self.centroid_x - BOUNDARY_LIMIT * self.x_max,
                self.centroid_x + BOUNDARY_LIMIT * self.x_max,
            ]
        )
        plt.ylim(
            [
                self.centroid_y - BOUNDARY_LIMIT * self.y_max,
                self.centroid_y + BOUNDARY_LIMIT * self.y_max,
            ]
        )

    def _on_press(self, event):
        """
        When the mouse is clicked the press variable gets assigned with current x and y coordinates.
        Input:
            event: click event information
        """
        self.press = event.xdata, event.ydata

    def _on_motion(self, event):
        """

        Input:
            event: Click event information
        """
        # This function gets called when the mouse moves
        if self.press is None:
            return  # return if mouse isnt pressed since the user cant drag anyway.
        xpress, ypress = self.press

        # Get mouse displacement by calculating the difference in current coordinates and previous coordinates.
        dx = event.xdata - xpress
        dy = event.ydata - ypress

        # Clear the current plotted data and get ready to plot new data.
        plt.cla()

        """
            Calculate how much we want to rotate along x and y axis separately.
            We calculate this by splitting 3 degrees(arbitrarly assigned valueof rotation per iteration)
            into 2 parts such that it is directly proportional to the projection of the displacement vector
            along the direction of movement of mouse.
        """
        if dx != 0.0 and dy != 0.0:
            angle = round((dx / (abs(dx) + abs(dy))) * TOTAL_ROTATION_PER_ITERATION, 3)
            self._rotate_y_3D(angle)
            angle = round((dy / (abs(dx) + abs(dy))) * TOTAL_ROTATION_PER_ITERATION, 3)
            self._rotate_x_3D(-angle)

        # reset the previous points variable
        self.press = event.xdata, event.ydata
        # Display newly calculated figure
        self.display()

    def _on_release(self, event):
        # On release of the mouse we reset the press data
        self.press = None

    def display(self):
        """
        Iterate over all the given faces of the figure. Each face will have 3 points (given constraint).
        Connect those 3 points to form a triangle, which will be one of the faces of our figure.
        """
        for face in self.faces:
            for i in range(len(face.coordinates)):
                for j in range(i + 1, len(face.coordinates)):
                    x_values = [
                        self.vertices[face.coordinates[i]].x_coordinate,
                        self.vertices[face.coordinates[j]].x_coordinate,
                    ]
                    y_values = [
                        self.vertices[face.coordinates[i]].y_coordinate,
                        self.vertices[face.coordinates[j]].y_coordinate,
                    ]
                    plt.plot(x_values, y_values, "bo", linestyle="-")

        # Set xlim, ylim of the plotted box so that figure doesnt go out of bounds.
        self._xy_boundary_limits()

        plt.pause(FIGURE_REFRESH_DELTA)

    def construct_figure(
        self,
        vertices: Dict[int, Tuple[float, float, float]],
        faces: List[Tuple[int, int, int]],
    ):
        """
        Initialize values.
        """

        # Step 1: Construct the vertices from the given coordinates
        for vertex_idx, point_coordinates in vertices.items():
            self.vertices[vertex_idx] = Point(*point_coordinates)

        # Step 2: Construct the faces using the vertices as a base and the
        # references of the points given in the input
        for face_coordinates in faces:
            self.faces.append(
                Face((face_coordinates[0], face_coordinates[1], face_coordinates[2]))
            )

        plt.rcParams["figure.figsize"] = [7.50, 7.50]
        plt.rcParams["figure.autolayout"] = True
        self._xy_max()
        self.display()
        plt.connect("button_press_event", self._on_press)
        plt.connect("button_release_event", self._on_release)
        plt.connect("motion_notify_event", self._on_motion)
        plt.show()
