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

    def __getitem__(self, x):
        return {0: self.x_coordinate, 1: self.y_coordinate, 2: self.z_coordinate}[x]


class Face:
    def __init__(
        self,
        coordinates: Tuple[Point, Point, Point],
        vertices: Dict[int, Tuple[float, float, float]],
    ):
        self.coordinates = coordinates
        self.vertices = vertices

    def get_centroid(self):

        cent_x = round(
            sum([self.vertices[idx].x_coordinate for idx in self.coordinates])
            / len(self.coordinates),
            3,
        )
        cent_y = round(
            sum([self.vertices[idx].y_coordinate for idx in self.coordinates])
            / len(self.coordinates),
            3,
        )
        cent_z = round(
            sum([self.vertices[idx].z_coordinate for idx in self.coordinates])
            / len(self.coordinates),
            3,
        )

        return [cent_x, cent_y, cent_z]

    def _get_point(self, idx):
        return self.vertices[self.coordinates[idx]]

    def get_face_plane_equation(self):
        p1 = self._get_point(0)
        p2 = self._get_point(1)
        p3 = self._get_point(2)

        a = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2])
        b = (p3[0] - p1[0]) * (p2[2] - p1[2]) - (p2[0] - p1[0]) * (p3[2] - p1[2])
        c = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])

        return a, b, c


class Figure:
    def __init__(self, num_vertices, num_faces):
        self.num_vertices = num_vertices
        self.num_faces = num_faces
        self.faces = []
        self.press = None
        self.centroid_x = None
        self.centroid_y = None
        self.centroid_z = None
        self.x_max = None
        self.y_max = None
        self.z_max = None
        self.vertices = {}

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
        z_coordinates = [vertex.z_coordinate for vertex in self.vertices.values()]

        self.centroid_x = sum(x_coordinates) / len(x_coordinates)
        self.centroid_y = sum(y_coordinates) / len(y_coordinates)
        self.centroid_z = sum(z_coordinates) / len(z_coordinates)

        self.x_max = max([abs(x) for x in x_coordinates])
        self.y_max = max([abs(y) for y in y_coordinates])
        self.z_max = max([abs(z) for z in z_coordinates])

    def _get_figure_centroids(self):
        x_coordinates = [vertex.x_coordinate for vertex in self.vertices.values()]
        y_coordinates = [vertex.y_coordinate for vertex in self.vertices.values()]
        z_coordinates = [vertex.z_coordinate for vertex in self.vertices.values()]

        return (
            sum(x_coordinates) / len(x_coordinates),
            sum(y_coordinates) / len(y_coordinates),
            sum(z_coordinates) / len(z_coordinates),
        )

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
        fig_cent_x, fig_cent_y, fig_cent_z = self._get_figure_centroids()

        for face in self.faces:

            # Get centroid of current face
            face_centroid = face.get_centroid()

            # Vector between centroids
            centroid_vector = [
                face_centroid[0] - fig_cent_x,
                face_centroid[1] - fig_cent_y,
                face_centroid[2] - fig_cent_z,
            ]

            # normal vectors to face plane
            a, b, c = face.get_face_plane_equation()
            plane_normal_1 = [a, b, c]
            plane_normal_2 = [-a, -b, -c]

            # Get dot product of both vectors of the face plane with centroid vector
            dot_product_pn1_centroid = sum(
                [t1 * t2 for t1, t2 in zip(centroid_vector, plane_normal_1)]
            )
            dot_product_pn2_centroid = sum(
                [t1 * t2 for t1, t2 in zip(centroid_vector, plane_normal_2)]
            )

            # Select the one which makes an acute angle with the centroid vector, this will be outwards facing.
            if dot_product_pn1_centroid > 0:
                outwards_facing_vector = plane_normal_1
            if dot_product_pn2_centroid >= 0:
                outwards_facing_vector = plane_normal_2

            # Calculate the cos angle value between outward facing normal vertor and the vector
            # directed towards viewer(+z axis, viewing direction)
            positive_z_axis_vector = [0, 0, 1]
            dot_product_normal_zaxis = sum(
                [
                    t1 * t2
                    for t1, t2 in zip(outwards_facing_vector, positive_z_axis_vector)
                ]
            ) / (
                math.sqrt(sum([t1 * t1 for t1 in outwards_facing_vector]))
                * math.sqrt(sum([t1 * t1 for t1 in positive_z_axis_vector]))
            )

            # Finally determine if face can be viewed by the user or not.
            if dot_product_normal_zaxis > 0:
                color_face = True
                color_of_face = dot_product_normal_zaxis  # since this is between 0 and 1 we can just use it as is.
            else:
                color_face = False

            if color_face:
                plt.fill(
                    "j",
                    "k",
                    color=(0.0, 0.0, color_of_face),
                    data={
                        "j": [
                            face._get_point(0)[0],
                            face._get_point(1)[0],
                            face._get_point(2)[0],
                        ],
                        "k": [
                            face._get_point(0)[1],
                            face._get_point(1)[1],
                            face._get_point(2)[1],
                        ],
                    },
                )

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
                Face(
                    (face_coordinates[0], face_coordinates[1], face_coordinates[2]),
                    self.vertices,
                )
            )

        plt.rcParams["figure.figsize"] = [7.50, 7.50]
        plt.rcParams["figure.autolayout"] = True
        self._xy_max()
        self.display()
        plt.connect("button_press_event", self._on_press)
        plt.connect("button_release_event", self._on_release)
        plt.connect("motion_notify_event", self._on_motion)
        plt.show()
