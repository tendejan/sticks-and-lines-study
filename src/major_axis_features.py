from scipy.spatial.transform import Rotation
import numpy as np

EULER_ORDER = "YZX" #as defined by the data 'coding' appilcation

VIEW_AXIS = np.array([0,0,-1]) #viewing from the z axis

def rotate_major_axis(axis_vector, euler_x, euler_y, euler_z):
    """rotate the main axis by the given euler sequence"""
    rotation = Rotation.from_euler(EULER_ORDER, (euler_y, euler_z, euler_x), degrees=True)
    rotated_axis = rotation.apply(axis_vector)
    return rotated_axis

def compute_projected_angle_from_vertical(major_axis, perceptual_upright=np.array([0,1])):
    """Projects the major axis to the xy plane and computes its angle from the perceptual upright"""
    projected_vector = major_axis[:2] # drop the z axis to project to xy plane
    #find the normalization factor
    norm_factor = np.linalg.norm(projected_vector)*np.linalg.norm(perceptual_upright)
    # take the dot product and divide by the norm factor
    shadow = np.dot(projected_vector, perceptual_upright)/norm_factor
    #find the angle from the absolute value of the normalized shadow
    angle_radians = np.acos(np.abs(shadow))
    return np.rad2deg(angle_radians) #convert to degrees and return
