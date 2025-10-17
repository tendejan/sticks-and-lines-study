import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Ellipse
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.lines import Line2D
from os.path import join
import polars as pl
from PIL import Image
import numpy as np

MAJOR_AXIS_RENDER_SCALE = 100

def render_props(dataframe, renditions_directory, visualization_directory):
    """save the rendered computed properties to a file and return an updated dataframe with the path to the rendered file"""
    visualized_image_paths = []
    fig, ax = plt.subplots(1, 1, figsize=(10,10))

    for row in dataframe.iter_rows(named=True):
        ax.clear()
        #extract the values from the dataframe
        rendition_path = join(renditions_directory, row['image'])
        x, y = row['centroid_x'], row['centroid_y']

        #this is the 3d projected vector
        dx = row['rotated_major_axis'][0] * MAJOR_AXIS_RENDER_SCALE
        dy = row['rotated_major_axis'][1] * -MAJOR_AXIS_RENDER_SCALE#needs to be fipped for a convention mismatch

        #these are the fit ellipse props
        ellipse_axis_major = row['axis_major_length']
        ellipse_axis_minor = row['axis_minor_length']
        ellipse_orientation = np.deg2rad(-row['orientation_degrees']) #needs to be flipped for a convention mismatch

        #display the renditon
        img = Image.open(rendition_path)
        ax.imshow(img)

        #draw the major axis
        arrow_patch = FancyArrow(x, y, dx, dy, color='r', width=7, head_width=15, head_length=10)
        ax.add_patch(arrow_patch)

        #draw the ellipse
        ellipse_patch = Ellipse(
            (x, y),
            2*ellipse_axis_major, 2*ellipse_axis_minor,
            angle=np.rad2deg(ellipse_orientation), #and this one uses degrees
            fill=False, edgecolor='b', linewidth=1,
        )
        ax.add_patch(ellipse_patch)

        #draw the ellipse moments
        # Calculate endpoints of unit-length major axis
        x1 = x - np.cos(ellipse_orientation) * ellipse_axis_major
        y1 = y - np.sin(ellipse_orientation) * ellipse_axis_major
        x2 = x + np.cos(ellipse_orientation) * ellipse_axis_major
        y2 = y + np.sin(ellipse_orientation) * ellipse_axis_major
        ellipse_axis_patch = Line2D(
            [x1,x2],
            [y1,y2],
            color='c', linewidth=3
        )
        ax.add_line(ellipse_axis_patch)

        #format the plot
        ax.set_title(row['image'])
        ax.axis('off')
        #Write the figure out
        out_path = join(visualization_directory, row['image'])
        visualized_image_paths.append(out_path)
        fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)

    #update the dataframe and return it
    dataframe_with_visuals = dataframe.with_columns(
        pl.Series('rendered_props_path', visualized_image_paths)
    )
    return dataframe_with_visuals

def visualize_fig1(combined_dataframe):
    # Create the figure with two subplots side by side
    fig, (ax_scatter, ax_image) = plt.subplots(1, 2, figsize=(20, 10), 
                                                gridspec_kw={'width_ratios': [1, 1]})

    # Create the scatter plot
    scatter = ax_scatter.scatter(
        combined_dataframe['projected_angle_from_vert'],
        combined_dataframe['orientation_degrees_from_vertical'],
        s=100,
        alpha=0.6
    )

    ax_scatter.set_xlabel('Angle of Major Axis Projection from Vertical')
    ax_scatter.set_ylabel('Orientation of Second Moment Ellipse')
    ax_scatter.set_title('Scatter Plot with Image Hover')

    # Set up the image axes
    ax_image.axis('off')
    ax_image.set_title('Hover over a point to see the image')
    
    # Initialize with empty image
    img_display = ax_image.imshow(np.ones((100, 100, 3)))
    
    # Text annotation for labels on scatter plot
    annot = ax_scatter.annotate(
        "",
        xy=(0, 0),
        xytext=(10, 10),
        textcoords="offset points",
        bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.9),
        fontsize=10,
        visible=False
    )

    def hover(event):
        if event.inaxes == ax_scatter:
            # Check if mouse is over a point
            cont, ind = scatter.contains(event)
            
            if cont:
                # Get the index of the point and convert to Python int
                idx = int(ind["ind"][0])
                
                # Get data for this point (Polars syntax)
                x = combined_dataframe['projected_angle_from_vert'][idx]
                y = combined_dataframe['orientation_degrees_from_vertical'][idx]
                image_name = combined_dataframe['image'][idx]
                image_path = combined_dataframe['rendered_props_path'][idx]
                
                # Update text annotation on scatter plot
                text = f"{image_name}\nMajor Axis: {x:.2f}\nEllipse: {y:.2f}"
                annot.xy = (x, y)
                annot.set_text(text)
                annot.set_visible(True)
                
                # Load and display image in right subplot
                try:
                    img = Image.open(image_path)
                    
                    # Update the image display
                    img_display.set_data(img)
                    img_display.set_extent([0, img.width, img.height, 0])
                    ax_image.set_xlim(0, img.width)
                    ax_image.set_ylim(img.height, 0)
                    ax_image.set_title(f'{image_name} - Major Axis: {x:.2f}° | Ellipse: {y:.2f}°')
                    
                except Exception as e:
                    print(f"Error loading image: {e}")
                    ax_image.set_title('Error loading image')
                
                fig.canvas.draw_idle()
            else:
                # Hide annotation when not hovering over points
                if annot.get_visible():
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # Connect the hover event
    fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.tight_layout()
    plt.show()