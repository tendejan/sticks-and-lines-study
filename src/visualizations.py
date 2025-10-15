import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Ellipse
from matplotlib.lines import Line2D
from PIL import Image
import numpy as np

def render_props(rendition_path, out_path, centroid_x, centriod_y, axis_major, axis_minor, ori_deg, rot_ax):
    """render the computed properties on top of the original image and save to a file"""
    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    MAJOR_AXIS_SCALE = 100

    #display the renditon
    img = Image.open(rendition_path)
    ax.imshow(img)

    x, y = centroid_x, centriod_y

    #draw the major axis
    dx = rot_ax[0] * MAJOR_AXIS_SCALE
    dy = rot_ax[1] * -MAJOR_AXIS_SCALE#needs to be fipped for a convention mismatch
    arrow_patch = FancyArrow(x, y, dx, dy, color='r', width=7, head_width=15, head_length=10)
    ax.add_patch(arrow_patch)

    #draw the ellipse
    ellipse_patch = Ellipse(
        (x, y),
        2*axis_major, 2*axis_minor,
        angle=ori_deg*-1, # needs to be flipped for convention mismatch
        fill=False, edgecolor='b', linewidth=1,
    )
    ax.add_patch(ellipse_patch)

    #draw the ellipse moments
    angle_rad = np.deg2rad(-ori_deg)
    # Calculate endpoints of unit-length major axis
    x1 = x - np.cos(angle_rad) * axis_major
    y1 = y - np.sin(angle_rad) * axis_major
    x2 = x + np.cos(angle_rad) * axis_major
    y2 = y + np.sin(angle_rad) * axis_major
    ellipse_axis_patch = Line2D(
        [x1,x2],
        [y1,y2],
        color='c', linewidth=3
    )
    ax.add_line(ellipse_axis_patch)

    #Write the figure out
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)

def visualize_fig1():
    ...