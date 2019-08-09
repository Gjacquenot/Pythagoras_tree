# -*- coding: UTF-8 -*-
# Works with Python2 and Python3
#
# Requires numpy and matplotlib to create pngs.
#
# Get started with : python Pythagoras_tree.py --help

def Pythagoras_tree(m=0.8, n=12):
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
    from math import pi, atan2, sqrt
    import numpy as np
    # Check inputs
    if m <= 0:
        raise Exception('Length of m has to be greater than zero')
    if int(n) != float(n):
        raise Exception('The number of level has to be integer')

    # Compute constants
    d = sqrt(1.0 + m ** 2)
    # Normalized length 1
    c1 = 1.0 / d
    # Normalized length 2
    c2 = m / d
    # Translation pattern
    T = np.array([[0.0, 1.0 / (1.0 + m ** 2)],
                  [1.0, 1.0 + m / (1.0 + m ** 2)]])
    # Defines the first rotation angle
    alpha1 = atan2(m, 1.0)
    # Defines the second rotation angle
    alpha2 = alpha1 - pi / 2.0
    # Defines pi2
    pi2 = 2.0 * pi
    # Number of elements (square)
    nEle = 2 ** (n + 1) - 1
    # Matrice containing the tree
    M = np.zeros((nEle, 5))
    # Initialization of the tree
    M[0, :] = [0.0, -1.0, 0.0, 1.0, 1.0]

    # Compute the level of each square contained in the resulting matrix
    offset = 0
    for i in range(n + 1):
        tmp = 2 ** i
        M[offset: offset + tmp, 4] = i
        offset += tmp

    def mat_rot(x):
        c = np.cos(x)
        s = np.sin(x)
        return np.array([[c, -s], [s, c]])
    # Compute the position and size of each square wrt its parent
    for i in range(1, nEle, 2):
        j = (i + 1) // 2 - 1
        mT = M[j, 3] * np.dot(mat_rot(M[j, 2]), T)
        Tx = mT[0, :] + M[j, 0]
        Ty = mT[1, :] + M[j, 1]
        theta1 = (M[j, 2] + alpha1) % pi2
        theta2 = (M[j, 2] + alpha2) % pi2
        M[i, 0: 4] = [Tx[0], Ty[0], theta1, M[j, 3] * c1]
        M[i + 1, 0: 4] = [Tx[1], Ty[1], theta2, M[j, 3] * c2]
    return M


def Pythagor_tree_write2svg(m=0.8, n=12, colormap='summer', M=[]):
    """
    A svg file giving a vectorial display of the tree. The name of
    file is generated from the parameter m, n, colormap. The file is
    stored in the current folder.
    """
    from math import pi, atan2, sqrt, ceil
    import numpy as np
    import matplotlib
    from matplotlib import cm
    Display_metadata = True

    nEle = M.shape[0]
    r2 = sqrt(2)
    LOffset = M[nEle - 1, 3] + 0.1
    min_x = np.min(M[:, 0] - r2 * M[:, 3]) - LOffset
    max_x = np.max(M[:, 0] + r2 * M[:, 3]) + LOffset
    min_y = np.min(M[:, 1]) - LOffset
    max_y = np.max(M[:, 1] + r2 * M[:, 3]) + LOffset

    # Compute the color of tree
    ColorM = cm.get_cmap(colormap)
    co = 100
    Wfig = ceil(co * (max_x - min_x))
    Hfig = ceil(co * (max_y - min_y))
    filename = 'Pythagoras_tree_' + str(m).replace('.', '_') + '__' + \
               str(n) + '__' + colormap + '.svg'
    fid = open(filename, 'wt')
    w = fid.write
    w('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    if not Display_metadata:
        w('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
        w('  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
    w('<svg width="{0}" height="{1}" version="1.1"\n'.format(Wfig, Hfig))
    # w(['<svg width="12cm" height="4cm" version="1.1"\n']); % Wfig,

    # w(['<svg width="15cm" height="10cm" '...
    #              'viewBox="0 0 %d %d" version="1.1"\n'],...
    #              Wfig,Hfig);
    if Display_metadata:
        w('\txmlns:dc="http://purl.org/dc/elements/1.1/"\n')
        w('\txmlns:cc="http://creativecommons.org/ns#"\n')
        w('\txmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
    w('\txmlns:svg="http://www.w3.org/2000/svg"\n')
    w('\txmlns="http://www.w3.org/2000/svg"\n')
    w('\txmlns:xlink="http://www.w3.org/1999/xlink">\n')

    if Display_metadata:
        w('\t<title>Pythagoras tree</title>\n')
        w('\t<metadata>\n')
        w('\t\t<rdf:RDF>\n')
        w('\t\t\t<cc:Work\n')
        w('\t\t\t\trdf:about="">\n')
        w('\t\t\t\t<dc:format>image/svg+xml</dc:format>\n')
        w('\t\t\t\t<dc:type\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://purl.org/dc/dcmitype/StillImage"/>\n')
        w('\t\t\t\t<dc:title>Pythagoras tree</dc:title>\n')
        w('\t\t\t\t<dc:creator>\n')
        w('\t\t\t\t\t<cc:Agent>\n')
        w('\t\t\t\t\t\t<dc:title>Guillaume Jacquenot</dc:title>\n')
        w('\t\t\t\t\t</cc:Agent>\n')
        w('\t\t\t\t</dc:creator>\n')
        w('\t\t\t\t<cc:license\n')
        w('\t\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/licenses/by-nc-sa/3.0/"/>\n')
        w('\t\t\t</cc:Work>\n')
        w('\t\t\t<cc:License\n')
        w('\t\t\t\trdf:about=' +
          '"http://creativecommons.org/licenses/by-nc-sa/3.0/">\n')
        w('\t\t\t\t<cc:permits\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#Reproduction"/>\n')
        w('\t\t\t\t<cc:permits\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#Reproduction"/>\n')
        w('\t\t\t\t<cc:permits\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#Distribution"/>\n')
        w('\t\t\t\t<cc:requires\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#Notice"/>\n')
        w('\t\t\t\t<cc:requires\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#Attribution"/>\n')
        w('\t\t\t\t<cc:prohibits\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#CommercialUse"/>\n')
        w('\t\t\t\t<cc:permits\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#DerivativeWorks"/>\n')
        w('\t\t\t\t<cc:requires\n')
        w('\t\t\t\t\trdf:resource=' +
          '"http://creativecommons.org/ns#ShareAlike"/>\n')
        w('\t\t\t</cc:License>\n')
        w('\t\t</rdf:RDF>\n')
        w('\t</metadata>\n')
    w('\t<defs>\n')
    w('\t\t<rect width="{0}" height="{1}" \n'.format(co, co))
    w('\t\t\tx="0" y="0"\n')
    w('\t\t\tstyle="fill-opacity:1;stroke:#00d900;stroke-opacity:1"\n')
    w('\t\t\tid="squa"\n')
    w('\t\t/>\n')
    w('\t</defs>\n')
    w('\t<g transform="translate({0} {1}) rotate(180) " >\n'.format(
        round(co * max_x), round(co * max_y)))
    for i in range(n + 1):
        w('\t\t<g style="fill:{0};" >\n'.format(
            matplotlib.colors.rgb2hex(ColorM(1.0 - float(i) / (n + 1)))))
        offset = 2**i - 1
        for j in range(2**i):
            k = j + offset
            w(('\t\t\t<use xlink:href="#squa" ' +
               'transform="translate({0:+010.5f} {1:+010.5f}) ' +
               'rotate({2:3.1f}) scale({3:8.6f})" />\n').format(
                  co * M[k, 0], co * M[k, 1], M[k, 2] * 180 / pi, M[k, 3]))
        w('\t\t</g>\n')
    w('\t</g>\n')
    w('</svg>\n')
    fid.close()


def Pythagor_tree_plot(M, colormap='summer', outputFilename='lm.png'):
    from math import pi
    from matplotlib import cm
    from matplotlib import patches
    import matplotlib.pyplot as plt
    ColorM = cm.get_cmap(colormap)
    fig, ax = plt.subplots()
    for i in range(M.shape[0]):
        cx = M[i, 0]
        cy = M[i, 1]
        theta = M[i, 2]
        si = M[i, 3]
        rect = patches.Rectangle([cx, cy], si, si,
                                 angle=theta * 180.0 / pi,
                                 ec="none",
                                 color=ColorM(1.0 - i / (M[-1, 4] + 1)))
        ax.add_patch(rect)
    plt.xlim([-4, 4])
    plt.ylim([-1.5, 3.5])
    # plt.gca().relim()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    fig.savefig(outputFilename, bbox_inches='tight')


def getListOfColormaps():
    LCmap = [
        'autumn',
        'bone',
        'colorcube',
        'cool',
        'copper',
        'flag',
        'gray',
        'hot',
        'hsv',
        'jet',
        'lines',
        'pink',
        'prism',
        'spring',
        'summer',
        'white',
        'winter']
    return LCmap


def isColormap(cmap):
    # This function returns true if 'cmap' is a valid colormap
    return cmap in getListOfColormaps()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='This script creates a SVG image of a Pythagoras tree, a' +
                    ' plane fractal constructed from squares.')
    pa = parser.add_argument
    pa('-r', '--ratio',
        type=float, default=1.0,
        help='r ( r > 0 ) is the relative length of one of the side ' +
             'right-angled triangle. ' +
             'The second side of the right-angle is taken to be one. ' +
             'To have a symmetric tree, r has to be 1.')
    pa('-l', '--level',
        type=int, default=10,
        help='n is the level of recursion. The number of elements of ' +
             'tree is equal to 2**(n+1)-1. A reasonable number for n ' +
             'is 10.')
    pa('-p', '--plot',
        action='store_true',
        help='Option used to display the tree as a PNG image with matplotlib')
    pa('-c', '--colormap',
        type=str, default='summer',
        help='Matplotlib colormap used to generate color of the different ' +
             'levels of the tree. ' +
             'Possible values are {0}'.format(' '.join(getListOfColormaps())))
    args = parser.parse_args()

    # Create a matrix containing the informations representing the tree
    # Each row represents a single square
    M = Pythagoras_tree(m=args.ratio, n=args.level)

    if args.plot:
        # Display the tree
        filename = 'Pythagoras_tree_' + \
            str(args.ratio).replace('.', '_') + '__' + \
            str(args.level) + '__' + args.colormap + '.png'
        Pythagor_tree_plot(M, colormap=args.colormap, outputFilename=filename)

    # Write results to an SVG file
    Pythagor_tree_write2svg(args.ratio, args.level, args.colormap, M)
