# -*- coding: UTF-8 -*-
# Works with Python2 and Python3

def Pythagoras_tree(m = 0.8, n = 12):
    # Compute Pythagoras_tree
    # The Pythagoras Tree is a plane fractal constructed from squares.
    # It is named after Pythagoras  because each triple of touching squares
    # encloses a right triangle, in a configuration traditionally used to
    # depict the Pythagorean theorem.
    # http://en.wikipedia.org/wiki/Pythagoras_tree
    #
    #       All these arguments are optional: the function can run with
    #       argument.
    # Output :
    #       - Matrix M: Pythagoras tree is stored in a matrix M.
    #         This matrix has 5 columns.
    #         Each row corresponds to the coordinate of each square of the tree
    #         The two first columns give the bottom-left position of each
    #         square. The third column corresponds to the orientation angle of
    #         each square. The fourth column gives the size of each square. The
    #         fifth column specifies the level of recursion of each square.
    #         The first row corresponds to the root of the tree. It is always
    #         the same
    #         M[0,:] = [0 -1 0 1 1];
    #         The leaf located at row i will give 2 leaves located at 2*i and
    #         2*i+1.

    from math import pi, atan2, sqrt
    import numpy as np
    # Check inputs
    if m <= 0:
        raise Exception('Length of m has to be greater than zero')
    if int(n)!=float(n):
        raise Exception('The number of level has to be integer')

    ## Compute constants
    d      = sqrt(1+m**2)                                   #
    c1     = 1/d                                            # Normalized length 1
    c2     = m/d                                            # Normalized length 2
    T      = np.array([[0,1/(1+m**2)],[1,1+m/(1+m**2)]])   # Translation pattern
    alpha1 = atan2(m,1)                                     # Defines the first rotation angle
    alpha2 = alpha1-pi/2                                    # Defines the second rotation angle
    pi2    = 2*pi                                           # Defines pi2
    nEle   = 2**(n+1)-1                                     # Number of elements (square)
    M      = np.zeros((nEle,5))                             # Matrice containing the tree
    M[0,:] = [0,-1,0,1,1]                                   # Initialization of the tree

    # Compute the level of each square contained in the resulting matrix
    Offset = 0
    for i in range(n+1):
        tmp = 2**i
        M[Offset : Offset + tmp, 4] = i
        Offset += tmp

    def mat_rot(x):
        c = np.cos(x)
        s = np.sin(x)
        return np.array([[c,-s],[s,c]])
    # Compute the position and size of each square wrt its parent
    for i in range(1,nEle,2):
        j          = (i+1)//2-1
        mT         = M[j,3] * np.dot(mat_rot(M[j,2]), T )
        Tx         = mT[0,:] + M[j,0]
        Ty         = mT[1,:] + M[j,1]
        theta1     = (M[j,2]+alpha1)%pi2
        theta2     = (M[j,2]+alpha2)%pi2
        M[i  ,0:4] = [Tx[0],Ty[0],theta1,M[j,3]*c1]
        M[i+1,0:4] = [Tx[1],Ty[1],theta2,M[j,3]*c2]
    return M

def Pythagor_tree_write2svg(m = 0.8, n = 12, colormap = 'summer', M = []):
    # A svg file giving a vectorial display of the tree. The name of
    # file is generated from the parameter m,n,colormap. The file is
    # stored in the current folder.
    #
    # Determine the bounding box of the tree with an offset
    # Display_metadata = false;
    from math import pi, atan2, sqrt, ceil
    import numpy as np
    import matplotlib
    from matplotlib import cm
    Display_metadata = True

    nEle    = M.shape[0]
    r2      = sqrt(2)
    LOffset = M[nEle-1,3] + 0.1
    min_x   = np.min(M[:,0]-r2*M[:,3]) - LOffset
    max_x   = np.max(M[:,0]+r2*M[:,3]) + LOffset
    min_y   = np.min(M[:,1]          ) - LOffset
    max_y   = np.max(M[:,1]+r2*M[:,3]) + LOffset

    # Compute the color of tree
    ColorM = cm.get_cmap(colormap)
    co   = 100;
    Wfig = ceil(co*(max_x-min_x))
    Hfig = ceil(co*(max_y-min_y))
    filename = 'Pythagoras_tree_'+str(m).replace('.','_')+'__'+str(n)+'__'+colormap+'.svg'
    fid  = open(filename, 'wt');
    fid.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    if not Display_metadata:
        fid.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
        fid.write('  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
    fid.write('<svg width="{0}" height="{1}" version="1.1"\n'.format(Wfig,Hfig))
    # fid.write(['<svg width="12cm" height="4cm" version="1.1"\n']); % Wfig,

    # fid.write(['<svg width="15cm" height="10cm" '...
    #              'viewBox="0 0 %d %d" version="1.1"\n'],...
    #              Wfig,Hfig);
    if Display_metadata:
        fid.write('\txmlns:dc="http://purl.org/dc/elements/1.1/"\n')
        fid.write('\txmlns:cc="http://creativecommons.org/ns#"\n')
        fid.write('\txmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
    fid.write('\txmlns:svg="http://www.w3.org/2000/svg"\n')
    fid.write('\txmlns="http://www.w3.org/2000/svg"\n')
    fid.write('\txmlns:xlink="http://www.w3.org/1999/xlink">\n')

    if Display_metadata:
        fid.write('\t<title>Pythagoras tree</title>\n')
        fid.write('\t<metadata>\n')
        fid.write('\t\t<rdf:RDF>\n')
        fid.write('\t\t\t<cc:Work\n')
        fid.write('\t\t\t\trdf:about="">\n')
        fid.write('\t\t\t\t<dc:format>image/svg+xml</dc:format>\n')
        fid.write('\t\t\t\t<dc:type\n')
        fid.write('\t\t\t\t\trdf:resource="http://purl.org/dc/dcmitype/StillImage" />\n')
        fid.write('\t\t\t\t<dc:title>Pythagoras tree</dc:title>\n')
        fid.write('\t\t\t\t<dc:creator>\n')
        fid.write('\t\t\t\t\t<cc:Agent>\n')
        fid.write('\t\t\t\t\t\t<dc:title>Guillaume Jacquenot</dc:title>\n')
        fid.write('\t\t\t\t\t</cc:Agent>\n')
        fid.write('\t\t\t\t</dc:creator>\n')
        fid.write('\t\t\t\t<cc:license\n')
        fid.write('\t\t\t\t\t\trdf:resource="http://creativecommons.org/licenses/by-nc-sa/3.0/" />\n')
        fid.write('\t\t\t</cc:Work>\n')
        fid.write('\t\t\t<cc:License\n')
        fid.write('\t\t\t\trdf:about="http://creativecommons.org/licenses/by-nc-sa/3.0/">\n')
        fid.write('\t\t\t\t<cc:permits\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#Reproduction" />\n')
        fid.write('\t\t\t\t<cc:permits\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#Reproduction" />\n')
        fid.write('\t\t\t\t<cc:permits\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#Distribution" />\n')
        fid.write('\t\t\t\t<cc:requires\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#Notice" />\n')
        fid.write('\t\t\t\t<cc:requires\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#Attribution" />\n')
        fid.write('\t\t\t\t<cc:prohibits\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#CommercialUse" />\n')
        fid.write('\t\t\t\t<cc:permits\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#DerivativeWorks" />\n')
        fid.write('\t\t\t\t<cc:requires\n')
        fid.write('\t\t\t\t\trdf:resource="http://creativecommons.org/ns#ShareAlike" />\n')
        fid.write('\t\t\t</cc:License>\n')
        fid.write('\t\t</rdf:RDF>\n')
        fid.write('\t</metadata>\n')
    fid.write('\t<defs>\n')
    fid.write('\t\t<rect width="{0}" height="{1}" \n'.format(co,co))
    fid.write('\t\t\tx="0" y="0"\n')
    fid.write('\t\t\tstyle="fill-opacity:1;stroke:#00d900;stroke-opacity:1"\n')
    fid.write('\t\t\tid="squa"\n')
    fid.write('\t\t/>\n')
    fid.write('\t</defs>\n')
    fid.write('\t<g transform="translate({0} {1}) rotate(180) " >\n'.format(round(co*max_x),round(co*max_y)))
    for i in range(n+1):
        fid.write('\t\t<g style="fill:{0};" >\n'.format(matplotlib.colors.rgb2hex(ColorM(1.0-i/(n+1)))))
        Offset = 2**i-1
        for j in range(2**i):
            k = j + Offset
            fid.write('\t\t\t<use xlink:href="#squa" transform="translate({0:+010.5f} {1:+010.5f}) rotate({2:3.1f}) scale({3:8.6f})" />\n'.format(co*M[k,0],co*M[k,1],M[k,2]*180/pi,M[k,3]))
        fid.write('\t\t</g>\n')
    fid.write('\t</g>\n')
    fid.write('</svg>\n')
    fid.close()

def Pythagor_tree_plot(M, colormap = 'summer', outputFilename = 'lm.png'):
    from math import pi
    from matplotlib import cm
    from matplotlib import patches
    import matplotlib.pyplot as plt
    ColorM = cm.get_cmap(colormap)
    fig, ax = plt.subplots()
    for i in range(M.shape[0]):
        cx    = M[i,0]
        cy    = M[i,1]
        theta = M[i,2]
        si    = M[i,3]
        rect = patches.Rectangle([cx,cy], si, si, angle = theta * 180.0/pi, ec="none", color = ColorM(1.0-i/(M[-1,4]+1)))
        ax.add_patch(rect)
    plt.xlim([-4, 4])
    plt.ylim([-1.5, 3.5])
    # plt.gca().relim()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    fig.savefig(outputFilename, bbox_inches='tight')

def getListOfColormaps():
    LCmap = [\
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
        description = 'This script creates a SVG image of a Pythagoras tree, a plane fractal constructed from squares.')
    parser.add_argument('-r','--ratio', type = float, default = 1.0,
                       help = 'r ( r > 0 ) is the relative length of one of the side right-angled triangle. '+
                              'The second side of the right-angle is taken to be one. '+
                              'To have a symmetric tree, r has to be 1.')
    parser.add_argument('-l','--level', type = int, default = 10,
                       help='n is the level of recursion. The number of elements of tree is equal '+
                            'to 2**(n+1)-1. A reasonable number for n is 10.')
    parser.add_argument('-p','--plot', action = 'store_true',
                       help='Option used to display the tree as a PNG image with matplotlib')
    parser.add_argument('-c','--colormap', type = str, default = 'summer',
                       help='Matplotlib colormap used to generate color of the different levels of the tree. '+
                       'Possible values are {0}'.format(' '.join(getListOfColormaps())))
    args = parser.parse_args()

    # Create a matrix containing the informations representing the tree
    # Each row represents a single square
    M = Pythagoras_tree(m = args.ratio, n = args.level)

    if args.plot:
        # Display the tree
        filename = 'Pythagoras_tree_'+\
            str(args.ratio).replace('.','_')+'__'+\
            str(args.level)+'__'+args.colormap+'.png'
        Pythagor_tree_plot(M, colormap = args.colormap, outputFilename = filename)

    # Write results to an SVG file
    Pythagor_tree_write2svg(args.ratio, args.level, args.colormap, M)

