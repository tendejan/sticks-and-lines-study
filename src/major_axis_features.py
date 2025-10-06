from scipy.spatial.transform import Rotation
import numpy as np

EULER_ORDER = "YZX"

VIEW_AXIS = np.array([0,0,-1]) #viewing from the z axis

def compute_main_axis_rotated(axis_vector, euler_x, euler_y, euler_z):
    """rotate the main axis by the given euler sequence"""
    rotation = Rotation.from_euler(EULER_ORDER, (euler_y, euler_z, euler_x), degrees=True)
    rotated_axis = rotation.apply(axis_vector)
    return rotated_axis

def compute_angle_from_vertical(major_axis, view_axis=VIEW_AXIS):
    projected_vector = major_axis - np.dot(major_axis, view_axis) * view_axis
    projected_vector_normalized = projected_vector / np.linalg.norm(projected_vector)
    angle_from_vert = np.arctan2(projected_vector_normalized[1], projected_vector_normalized[0])
    angle_from_vert_deg = np.rad2deg(angle_from_vert)
    angle_from_vert_deg = angle_from_vert_deg % 90
    angle_from_vert_deg = min(angle_from_vert_deg, 90 + angle_from_vert_deg)
    return angle_from_vert_deg

def compute_angle_from_vertical_dot(major_axis, perceptual_upright=np.array([0,1])):
    projected_vector = major_axis[:2] # drop the z axis to project to xy plane
    #take the dot product to get the angle
    angle_radians = np.acos((np.abs(np.dot(projected_vector, perceptual_upright)))/(np.linalg.norm(projected_vector)*np.linalg.norm(perceptual_upright)))
    return np.rad2deg(angle_radians)