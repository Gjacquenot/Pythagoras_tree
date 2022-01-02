"""
This Python3 script creates a Pythagoras tree in SVG or PNG format

Requires numpy and matplotlib to create pngs.

Get started with : python Pythagoras_tree.py --help
"""
import argparse
from math import atan2, ceil, pi, sqrt
from typing import List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, patches


def pythagoras_tree(ratio: float = 1.0, nb_levels: int = 12):
    """
    Compute Pythagoras_tree
    The Pythagoras Tree is a plane fractal constructed from squares.
    It is named after Pythagoras because each triple of touching squares
    encloses a right triangle, in a configuration traditionally used to
    depict the Pythagorean theorem.
    http://en.wikipedia.org/wiki/Pythagoras_tree

          All these arguments are optional: the function can run with
          argument.
    Output :
          - Matrix M: Pythagoras tree is stored in a matrix M.
            This matrix has 5 columns.
            Each row corresponds to the coordinate of each square of the tree
            The two first columns give the bottom-left position of each
            square. The third column corresponds to the orientation angle of
            each square. The fourth column gives the size of each square. The
            fifth column specifies the level of recursion of each square.
            The first row corresponds to the root of the tree. It is always
            the same
            M[0,:] = [0 -1 0 1 1];
            The leaf located at row i will give 2 leaves located at 2*i and
            2*i+1.
    """
    # pylint: disable=too-many-locals
    # Check inputs
    if ratio <= 0:
        raise Exception("Length of ratio has to be greater than zero")
    if int(nb_levels) != float(nb_levels):
        raise Exception("The number of level has to be integer")

    # Compute constants
    c_d = sqrt(1.0 + ratio ** 2)
    # Normalized length 1
    c_1 = 1.0 / c_d
    # Normalized length 2
    c_2 = ratio / c_d
    # Translation pattern
    tr_pat = np.array(
        [[0.0, 1.0 / (1.0 + ratio ** 2)], [1.0, 1.0 + ratio / (1.0 + ratio ** 2)]]
    )
    # Defines the first rotation angle
    alpha1 = atan2(ratio, 1.0)
    # Defines the second rotation angle
    alpha2 = alpha1 - pi / 2.0
    # Number of elements (square)
    nb_elements = 2 ** (nb_levels + 1) - 1
    # Matrice containing the tree
    pt_arrray = np.zeros((nb_elements, 5))
    # Initialization of the tree
    pt_arrray[0, :] = [0.0, -1.0, 0.0, 1.0, 1.0]

    # Compute the level of each square contained in the resulting matrix
    offset = 0
    for i in range(nb_levels + 1):
        tmp = 2 ** i
        pt_arrray[offset : offset + tmp, 4] = i
        offset += tmp

    def mat_rot(angle_rad: float) -> np.ndarray:
        c_a = np.cos(angle_rad)
        s_a = np.sin(angle_rad)
        return np.array([[c_a, -s_a], [s_a, c_a]])

    # Compute the position and size of each square wrt its parent
    for i in range(1, nb_elements, 2):
        j = (i + 1) // 2 - 1
        t_m = pt_arrray[j, 3] * mat_rot(pt_arrray[j, 2]) @ tr_pat
        t_x = t_m[0, :] + pt_arrray[j, 0]
        t_y = t_m[1, :] + pt_arrray[j, 1]
        theta1 = (pt_arrray[j, 2] + alpha1) % (2.0 * pi)
        theta2 = (pt_arrray[j, 2] + alpha2) % (2.0 * pi)
        pt_arrray[i, 0:4] = [t_x[0], t_y[0], theta1, pt_arrray[j, 3] * c_1]
        pt_arrray[i + 1, 0:4] = [t_x[1], t_y[1], theta2, pt_arrray[j, 3] * c_2]
    return pt_arrray


def _svg_write_metadata(write):
    write("\t<title>Pythagoras tree</title>\n")
    write("\t<metadata>\n")
    write("\t\t<rdf:RDF>\n")
    write("\t\t\t<cc:Work\n")
    write('\t\t\t\trdf:about="">\n')
    write("\t\t\t\t<dc:format>image/svg+xml</dc:format>\n")
    write("\t\t\t\t<dc:type\n")
    write("\t\t\t\t\trdf:resource=" + '"http://purl.org/dc/dcmitype/StillImage"/>\n')
    write("\t\t\t\t<dc:title>Pythagoras tree</dc:title>\n")
    write("\t\t\t\t<dc:creator>\n")
    write("\t\t\t\t\t<cc:Agent>\n")
    write("\t\t\t\t\t\t<dc:title>Guillaume Jacquenot</dc:title>\n")
    write("\t\t\t\t\t</cc:Agent>\n")
    write("\t\t\t\t</dc:creator>\n")
    write("\t\t\t\t<cc:license\n")
    write(
        "\t\t\t\t\t\trdf:resource="
        + '"http://creativecommons.org/licenses/by-nc-sa/3.0/"/>\n'
    )
    write("\t\t\t</cc:Work>\n")
    write("\t\t\t<cc:License\n")
    write(
        "\t\t\t\trdf:about=" + '"http://creativecommons.org/licenses/by-nc-sa/3.0/">\n'
    )
    write("\t\t\t\t<cc:permits\n")
    write(
        "\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#Reproduction"/>\n'
    )
    write("\t\t\t\t<cc:permits\n")
    write(
        "\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#Reproduction"/>\n'
    )
    write("\t\t\t\t<cc:permits\n")
    write(
        "\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#Distribution"/>\n'
    )
    write("\t\t\t\t<cc:requires\n")
    write("\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#Notice"/>\n')
    write("\t\t\t\t<cc:requires\n")
    write("\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#Attribution"/>\n')
    write("\t\t\t\t<cc:prohibits\n")
    write(
        "\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#CommercialUse"/>\n'
    )
    write("\t\t\t\t<cc:permits\n")
    write(
        "\t\t\t\t\trdf:resource="
        + '"http://creativecommons.org/ns#DerivativeWorks"/>\n'
    )
    write("\t\t\t\t<cc:requires\n")
    write("\t\t\t\t\trdf:resource=" + '"http://creativecommons.org/ns#ShareAlike"/>\n')
    write("\t\t\t</cc:License>\n")
    write("\t\t</rdf:RDF>\n")
    write("\t</metadata>\n")


def pythagor_tree_write2svg(
    pt_array: np.ndarray,
    ratio: float = 0.8,
    nb_levels: int = 12,
    colormap_name: str = "summer",
):
    """
    A svg file giving a vectorial display of the tree. The name of
    file is generated from the parameter m, n, colormap. The file is
    stored in the current folder.
    """
    # pylint: disable=too-many-locals
    display_metadata = True

    nb_elements = pt_array.shape[0]
    length_offset = pt_array[nb_elements - 1, 3] + 0.1
    min_x = np.min(pt_array[:, 0] - sqrt(2) * pt_array[:, 3]) - length_offset
    max_x = np.max(pt_array[:, 0] + sqrt(2) * pt_array[:, 3]) + length_offset
    min_y = np.min(pt_array[:, 1]) - length_offset
    max_y = np.max(pt_array[:, 1] + sqrt(2) * pt_array[:, 3]) + length_offset

    # Compute the color of tree
    colormap = cm.get_cmap(colormap_name)
    nb_pixels = 100
    fig_w = ceil(nb_pixels * (max_x - min_x))
    fig_h = ceil(nb_pixels * (max_y - min_y))
    filename = (
        "Pythagoras_tree_"
        + str(ratio).replace(".", "_")
        + "__"
        + str(nb_levels)
        + "__"
        + colormap_name
        + ".svg"
    )
    with open(filename, "wt", encoding="utf-8") as fid:
        write = fid.write
        write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
        if not display_metadata:
            write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
            write('  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
        write('<svg width="{0}" height="{1}" version="1.1"\n'.format(fig_w, fig_h))
        # w(['<svg width="12cm" height="4cm" version="1.1"\n']); % Wfig,

        # w(['<svg width="15cm" height="10cm" '...
        #              'viewBox="0 0 %d %d" version="1.1"\n'],...
        #              Wfig,Hfig);
        if display_metadata:
            write('\txmlns:dc="http://purl.org/dc/elements/1.1/"\n')
            write('\txmlns:cc="http://creativecommons.org/ns#"\n')
            write('\txmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
        write('\txmlns:svg="http://www.w3.org/2000/svg"\n')
        write('\txmlns="http://www.w3.org/2000/svg"\n')
        write('\txmlns:xlink="http://www.w3.org/1999/xlink">\n')

        if display_metadata:
            _svg_write_metadata(write)
        write("\t<defs>\n")
        write('\t\t<rect width="{0}" height="{0}" \n'.format(nb_pixels))
        write('\t\t\tx="0" y="0"\n')
        write('\t\t\tstyle="fill-opacity:1;stroke:#00d900;stroke-opacity:1"\n')
        write('\t\t\tid="squa"\n')
        write("\t\t/>\n")
        write("\t</defs>\n")
        write(
            '\t<g transform="translate({0} {1}) rotate(180) " >\n'.format(
                round(nb_pixels * max_x), round(nb_pixels * max_y)
            )
        )
        for i in range(nb_levels + 1):
            write(
                '\t\t<g style="fill:{0};" >\n'.format(
                    matplotlib.colors.rgb2hex(
                        colormap(1.0 - float(i) / (nb_levels + 1))
                    )
                )
            )
            offset = 2 ** i - 1
            for j in range(2 ** i):
                k = j + offset
                write(
                    (
                        '\t\t\t<use xlink:href="#squa" '
                        + 'transform="translate({0:+010.5f} {1:+010.5f}) '
                        + 'rotate({2:3.1f}) scale({3:8.6f})" />\n'
                    ).format(
                        nb_pixels * pt_array[k, 0],
                        nb_pixels * pt_array[k, 1],
                        pt_array[k, 2] * 180 / pi,
                        pt_array[k, 3],
                    )
                )
            write("\t\t</g>\n")
        write("\t</g>\n")
        write("</svg>\n")


def pythagor_tree_plot(
    pt_array: np.ndarray, colormap_name: str = "summer", output_filename: str = "lm.png"
):
    """Plot a Pythagoras tree for a PNG format"""
    colormap = cm.get_cmap(colormap_name)
    fig, axis = plt.subplots()
    for i in range(pt_array.shape[0]):
        c_x = pt_array[i, 0]
        c_y = pt_array[i, 1]
        theta = pt_array[i, 2] * 180.0 / pi
        s_i = pt_array[i, 3]
        rect = patches.Rectangle(
            [c_x, c_y],
            s_i,
            s_i,
            angle=theta,
            ec="none",
            color=colormap(1.0 - i / (pt_array[-1, 4] + 1)),
        )
        axis.add_patch(rect)
    plt.xlim([-4, 4])
    plt.ylim([-1.5, 3.5])
    # plt.gca().relim()
    plt.gca().set_aspect("equal", adjustable="box")
    plt.axis("off")
    fig.savefig(output_filename, bbox_inches="tight", dpi=300)


def get_list_of_colormaps() -> List[str]:
    """Return the list of avalaible colors"""
    return [
        "autumn",
        "bone",
        "colorcube",
        "cool",
        "copper",
        "flag",
        "gray",
        "hot",
        "hsv",
        "jet",
        "lines",
        "pink",
        "prism",
        "spring",
        "summer",
        "white",
        "winter",
    ]


def is_colormap(cmap: str) -> bool:
    """This function returns True if 'cmap' is a valid colormap"""
    return cmap in get_list_of_colormaps()


def main():
    """Main entrypoint"""
    parser = argparse.ArgumentParser(
        description="This script creates a SVG image of a Pythagoras tree, a"
        + " plane fractal constructed from squares."
    )
    parser.add_argument(
        "-r",
        "--ratio",
        type=float,
        default=1.0,
        help="r ( r > 0 ) is the relative length of one of the side "
        + "right-angled triangle. "
        + "The second side of the right-angle is taken to be one. "
        + "To have a symmetric tree, r has to be 1.",
    )
    parser.add_argument(
        "-l",
        "--level",
        type=int,
        default=10,
        help="n is the level of recursion. The number of elements of "
        + "tree is equal to 2**(n+1)-1. A reasonable number for n "
        + "is 10.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="Option used to display the tree as a PNG image with matplotlib",
    )
    parser.add_argument(
        "-c",
        "--colormap",
        type=str,
        default="summer",
        help="Matplotlib colormap used to generate color of the different "
        + "levels of the tree. "
        + "Possible values are {0}".format(" ".join(get_list_of_colormaps())),
    )
    args = parser.parse_args()

    # Create a matrix containing the informations representing the tree
    # Each row represents a single square
    pt_array = pythagoras_tree(ratio=args.ratio, nb_levels=args.level)

    if args.plot:
        # Display the tree
        filename = (
            "Pythagoras_tree_"
            + str(args.ratio).replace(".", "_")
            + "__"
            + str(args.level)
            + "__"
            + args.colormap
            + ".png"
        )
        pythagor_tree_plot(
            pt_array, colormap_name=args.colormap, output_filename=filename
        )

    # Write results to an SVG file
    pythagor_tree_write2svg(pt_array, args.ratio, args.level, args.colormap)


if __name__ == "__main__":
    main()
