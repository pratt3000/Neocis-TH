import math
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from typing import Dict, List, AnyStr
from utils import driver

vertices = dict()
faces = list()
num_faces = 0
num_vertices = 0
xmin, xmax, ymin, ymax = 0, 0, 0, 0


def rotate_x_3D(theta: float, vertices: Dict[int, List[float]]):
    """
    Rotates all of the vertices of the given figure according to the theta given along the 'x-axis'
        Input:
            theta: The angle by which we need to rotate
            vertices: The coordinates of all vertices of the figure
        Output:
            vertices: updated vertices after rotating the coordinate points by theta
    """
    sinTheta = math.sin(theta)
    cosTheta = math.cos(theta)

    # Rotating all of the vectors of depending on the theta given using basic projectio formula.
    for point in vertices:
        x, y, z = vertices[point]

        cur_x = x
        cur_y = y * cosTheta - z * sinTheta
        cur_z = z * cosTheta + y * sinTheta

        vertices[point] = [cur_x, cur_y, cur_z]

    return vertices


def rotate_y_3D(theta: float, vertices: Dict[int, List[float]]):
    """
    Rotates all of the vertices of the given figure according to the theta given along the 'y-axis'
        Input:
            theta: The angle by which we need to rotate
            vertices: The coordinates of all vertices of the figure
        Output:
            vertices: updated vertices after rotating the coordinate points by theta
    """
    sinTheta = math.sin(theta)
    cosTheta = math.cos(theta)

    # Rotating all of the vectors of depending on the theta given using basic projectio formula.
    for point in vertices:
        x, y, z = vertices[point]

        cur_x = x * cosTheta + z * sinTheta
        cur_y = y
        cur_z = z * cosTheta - x * sinTheta

        vertices[point] = [cur_x, cur_y, cur_z]

    return vertices


def rotate_z_3D(theta: float, vertices: Dict[int, List[float]]):
    """
    *** We dont really need/use this function but this is there just for symmetry and incase its needed in the future ***
    Rotates all of the vertices of the given figure according to the theta given along the 'z-axis'
        Input:
            theta: The angle by which we need to rotate
            vertices: The coordinates of all vertices of the figure
        Output:
            vertices: updated vertices after rotating the coordinate points by theta
    """
    sinTheta = math.sin(theta)
    cosTheta = math.cos(theta)

    # Rotating all of the vectors of depending on the theta given using basic projectio formula.
    for point in vertices:
        x, y, z = vertices[point]

        cur_x = x * cosTheta - y * sinTheta
        cur_y = y * cosTheta + x * sinTheta
        cur_z = z

        vertices[point] = [cur_x, cur_y, cur_z]

    return vertices


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

    return num_vertices, num_faces, vertices, faces


def get_max_xy(vertices):
    """
    This function returns the xlim and ylim of the plot box so that the diagram doesnt go out of bounds.
    The coordinates of centroid of all points are calculated and then we leave space of 1.5 times the
    absoluate max_value that we found along that particular axis just to ensure that when the diagram
    rotates it can never go out of bounds.
        Input:
            vertices: The coordinates of all vertices of the figure.
            faces: The points which apppear on each face are given.
        Output:
            x_lowerbound: lower bound for x value of box-plot.
            x_upperbound: upper bound for x value of box-plot.
            y_lowerbound: lower bound for y value of box-plot.
            y_upperbound: upper bound for y value of box-plot.
    """
    # These lists maintain the x and y coordinates of all points.
    all_x = []
    all_y = []

    for point in vertices:
        all_x.append(vertices[point][0])
        all_y.append(vertices[point][1])

    # Calculating centroid of all points in figure.
    centroidx = sum(all_x) / len(all_x)
    centroidy = sum(all_y) / len(all_y)

    # Max of all absolute values along an axis.
    max_x = max([abs(x) for x in all_x])
    max_y = max([abs(y) for y in all_y])

    # Get lower and upper bounds of the box-plot.
    x_lowerbound = centroidx - 1.5 * max_x
    x_upperbound = centroidx + 1.5 * max_x
    y_lowerbound = centroidy - 1.5 * max_y
    y_upperbound = centroidy + 1.5 * max_y

    return x_lowerbound, x_upperbound, y_lowerbound, y_upperbound


def get_centroid_of_figure(vertices):
    # These lists maintain the x and y coordinates of all points.
    all_x = []
    all_y = []
    all_z = []

    for point in vertices:
        all_x.append(vertices[point][0])
        all_y.append(vertices[point][1])
        all_z.append(vertices[point][2])

    # Calculating centroid of all points in figure.
    centroid_x = round(sum(all_x) / len(all_x), 3)
    centroid_y = round(sum(all_y) / len(all_y), 3)
    centroid_z = round(sum(all_z) / len(all_z), 3)

    return [centroid_x, centroid_y, centroid_z]


def get_centroid_of_face(p1, p2, p3):

    cent_x = round((p1[0] + p2[0] + p3[0]) / 3, 3)
    cent_y = round((p1[1] + p2[1] + p3[1]) / 3, 3)
    cent_z = round((p1[2] + p2[2] + p3[2]) / 3, 3)
    return [cent_x, cent_y, cent_z]


def get_face_plane_equation(p1, p2, p3):
    # plane eqn: ax+by+cz+d=0
    a = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2])
    b = (p3[0] - p1[0]) * (p2[2] - p1[2]) - (p2[0] - p1[0]) * (p3[2] - p1[2])
    c = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])

    return a, b, c


def display_current_plot(vertices, faces):
    """
    Display the plot given the coordinates of current set of vertices.
    This function is called iterativey as we rotate the figure using GUI.
        Input:
            vertices: The coordinates of all vertices of the figure.
            faces: The points which apppear on each face are given.
        Output: None
    """
    fig_centroid = get_centroid_of_figure(vertices)

    # Go over all of the faces given in the object file and plot lines.
    for points in faces:

        # Given Constraint: Only 3 points per plane are present.
        # Get the (x,y,z coordinates of those points)
        p1 = [vertices[points[0]][0], vertices[points[0]][1], vertices[points[0]][2]]
        p2 = [vertices[points[1]][0], vertices[points[1]][1], vertices[points[1]][2]]
        p3 = [vertices[points[2]][0], vertices[points[2]][1], vertices[points[2]][2]]

        # Get centroid of current face
        face_centroid = get_centroid_of_face(p1, p2, p3)

        # vector between centroids
        centroid_vector = [t1 - t2 for t1, t2 in zip(face_centroid, fig_centroid)]

        # normal vectors to face plane
        a, b, c = get_face_plane_equation(p1, p2, p3)
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

        # Calculate the cos angle value between outward faceing normal vertor and vector facing viewer
        positive_z_axis_vector = [0, 0, 1]
        dot_product_normal_zaxis = sum(
            [t1 * t2 for t1, t2 in zip(outwards_facing_vector, positive_z_axis_vector)]
        ) / (
            math.sqrt(sum([t1 * t1 for t1 in outwards_facing_vector]))
            * math.sqrt(sum([t1 * t1 for t1 in positive_z_axis_vector]))
        )

        if dot_product_normal_zaxis > 0:
            color_face = True
            color_of_face = dot_product_normal_zaxis
        else:
            color_face = False

        if color_face:
            plt.fill(
                "j",
                "k",
                color=(0.0, 0.0, color_of_face),
                data={"j": [p1[0], p2[0], p3[0]], "k": [p1[1], p2[1], p3[1]]},
            )

    # Set xlim, ylim of the plotted box so that figure doesnt go out of bounds.
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.pause(0.5)


class DraggableFigure:
    def __init__(self, faces, vertices):
        self.faces = faces
        self.vertices = vertices
        self.press = None

    def connect(self):
        # Connect to all the events we need
        plt.connect("button_press_event", self.on_press)
        plt.connect("button_release_event", self.on_release)
        plt.connect("motion_notify_event", self.on_motion)

    def on_press(self, event):
        # Set press to datapoints since the mouse is being pressed and thus we may rotate the figure.
        self.press = event.xdata, event.ydata

    def on_motion(self, event):
        global vertices
        # This function gets called when the mouse moves
        if self.press is None:
            return  # return if mouse isnt pressed since the user isnt dragging.
        xpress, ypress = self.press

        # Get mouse displacement by calculating the difference in current coordinates and previous coordinates.
        dx = event.xdata - xpress
        dy = event.ydata - ypress

        # clear the current plotted data and get ready to plot new data.
        plt.cla()

        """ 
        Calculate how much we want to rotate along x and y axis separately.
        We calculate this by splitting 3 degrees(arbitrarly assigned valueof rotation per iteration)
        into 2 parts such that it is directly proportional to the projection of the displacement vector
        along the direction of movement of mouse.
        """
        if dx != 0.0 and dy != 0.0:
            angle = round((dx / (abs(dx) + abs(dy))) * 0.32, 3)
            vertices = rotate_y_3D(angle, vertices)
            angle = round((dy / (abs(dx) + abs(dy))) * 0.32, 3)
            vertices = rotate_x_3D(-angle, vertices)

        # reset the previous points variable
        self.press = event.xdata, event.ydata

        # Display newly calculated figure
        display_current_plot(vertices, self.faces)

    def on_release(self, event):
        # On release of the mouse we reset the press data
        self.press = None


def main():
    global vertices
    global faces
    global num_vertices
    global num_faces
    global xmin
    global xmax
    global ymin
    global ymax

    # Read file content
    num_vertices, num_faces, vertices, faces = read_object_file(path="object.txt")
    print(vertices)

    plt.rcParams["figure.figsize"] = [7.50, 7.50]
    plt.rcParams["figure.autolayout"] = True
    xmin, xmax, ymin, ymax = get_max_xy(vertices)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])

    display_current_plot(vertices, faces)

    dr = DraggableFigure(faces, vertices)
    dr.connect()
    plt.show()


if __name__ == "__main__":
    # main()
    driver()
