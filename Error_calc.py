import numpy as np
import geopandas as gpd
import scipy.io
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from shapely.geometry import Point, Polygon  # Import Polygon here
import matplotlib.colors as mcolors

# Define a custom colormap from green to white to red
cmap = mcolors.LinearSegmentedColormap.from_list(
    'green_white_red', [(0, 'green'), (0.5, 'orange'), (1, 'red')]
)
# Load data
file_path = 'synthetic_model2_T_3_res0.3.mat'
recover_file = 'recovered_model2_T_3_res0.3.mat'
cell_size_km = 0.2
area_analisis_km_square = 472
gdf = gpd.read_file("Natural_park.kml")  # Load KML with the park boundary
min_lon = 99.1319999999999908
min_lat = 1.4640000000000002
max_lon = 99.3389999999999986
max_lat = 1.6759999999999999
# Load station data
sta = np.loadtxt('stations.csv', usecols=0, skiprows=1, dtype=str, delimiter=';')
sta_lat = np.loadtxt('stations.csv', usecols=2, skiprows=1, delimiter=';')
sta_lon = np.loadtxt('stations.csv', usecols=1, skiprows=1, delimiter=';')

# Load the .mat files
mat_data = scipy.io.loadmat(file_path)
Synth = mat_data['V_2D']
grid_x = mat_data['x_grid_eff'][0]
grid_y = mat_data['y_grid_eff'][0]

grid_x = np.linspace(min_lon, max_lon, grid_x.shape[0])
grid_y = np.linspace(min_lat, max_lat, grid_y.shape[0])

mat_data = scipy.io.loadmat(recover_file)
Recovered = mat_data['V_map']


# Create the convex hull from the station locations (stations network boundary)
station_coords = np.column_stack((sta_lon, sta_lat))  # Combine lat and lon into a single array
hull = ConvexHull(station_coords)

# Get the polygon vertices (the boundary of the convex hull)
hull_vertices = station_coords[hull.vertices]

# Check if a point is inside the convex hull using the convex hull method
def is_point_inside_convex_hull(lon, lat, hull_vertices):
    point = Point(lon, lat)
    # Check if the point is inside the convex hull
    return point.within(Polygon(hull_vertices))

# Calculate the area of the convex hull
hull_polygon = Polygon(hull_vertices)
area_covered_by_hull = hull_polygon.area

# Print the area (in square degrees or in a different unit depending on your coordinate system)
print(f"Area covered by the convex hull: {area_covered_by_hull} square degrees")

# Loop through each grid point and check if it is inside the network (convex hull)
for i in range(Recovered.shape[0]):  # Loop over rows (latitudes)
    for j in range(Recovered.shape[1]):  # Loop over columns (longitudes)
        # If grid_x and grid_y are 1D, we just access them as grid_x[j] and grid_y[i]
        lon = grid_x[j] if grid_x.ndim == 1 else grid_x[0, j]
        lat = grid_y[i] if grid_y.ndim == 1 else grid_y[i, 0]

        # Check if the grid point is inside the convex hull
        if not is_point_inside_convex_hull(lon, lat, hull_vertices):
            # Assign NaN to points outside the network
            continue
            #Recovered[i, j] = np.nan
            #Synth[i, j] = np.nan

#---
# Plotting
fig, ax = plt.subplots(ncols=2, figsize=(10, 8))

# Plot the error using a pcolormesh
em = ax[0].pcolormesh(np.linspace(np.min(min_lon), np.max(max_lon), np.shape(grid_x)[0]),
                   np.linspace(np.min(min_lat), np.max(max_lat), np.shape(grid_y)[0]),
                   Synth, shading='auto', cmap='bwr')

em = ax[1].pcolormesh(np.linspace(np.min(min_lon), np.max(max_lon), np.shape(grid_x)[0]),
                   np.linspace(np.min(min_lat), np.max(max_lat), np.shape(grid_y)[0]),
                   Recovered.T, shading='auto', cmap='bwr', vmin=np.nanmin(Synth), vmax=np.nanmax(Synth))

# Plot the boundary of the Natural Park
gdf.plot(ax=ax[0], facecolor='none', edgecolor='black', linewidth=2)
gdf.plot(ax=ax[1], facecolor='none', edgecolor='black', linewidth=2)

# Plot the stations
ax[0].plot(sta_lon, sta_lat, 'k^', mec='white')
ax[1].plot(sta_lon, sta_lat, 'k^', mec='white')

# Customize the colorbar
#cbar = plt.colorbar(em, ax=ax[1])

# Set custom ticks and labels for the colorbar
#cbar.set_ticks([0, 25, 50])  # Ticks at Low (0), Medium, and High
#cbar.set_ticklabels(['High (Error 0 %)', 'Medium (Error 25 %)', 'Low (Error 50 %)'])
# Move the colorbar label to the top
#cbar.ax[0].set_title("Confidence")
# Labels and title
ax[0].set_xlabel("Longitude (º)")
ax[1].set_xlabel("Longitude (º)")
ax[0].set_ylabel("Latitude (º)")

ax[0].set_title('Synthetic')
ax[1].set_title('Recovered')

# Save the figure
plt.savefig('Synth_Recovered.png', dpi=600)
plt.savefig('Synth_Recovered.pdf', dpi=600, format='PDF')

# Show the plot
plt.show()

# Calculate the error as the relative error between recovered and synthetic models
Recovered = Recovered.T
Synth = Synth.T
Error = (Synth - Recovered)/Synth
#Error[Error <= -0.01] = Error*-1
Error = np.abs(Error) * 100

# Define a threshold for "close to 0" error (e.g., 1% error)
threshold = 10.0  # You can adjust this as needed

# Count the number of pixels with error less than or equal to the threshold
close_to_zero_count = np.sum(Error <= threshold)

# Calculate the total number of pixels
total_pixels = Error.size

# Calculate the percentage of pixels with error close to 0
percentage_close_to_zero = (close_to_zero_count / total_pixels) * 100

print('Porcentage with high confidence', percentage_close_to_zero)
print(close_to_zero_count)
print(area_covered_by_hull)
print('Area covered', (percentage_close_to_zero * area_analisis_km_square)/100)


# Plotting
fig, ax = plt.subplots(figsize=(10, 8))

# Plot the error using a pcolormesh

em = ax.pcolormesh(np.linspace(np.min(min_lon), np.max(max_lon), np.shape(grid_x)[0]),
                   np.linspace(np.min(min_lat), np.max(max_lat), np.shape(grid_y)[0]),
                   Error, shading='auto', cmap=cmap,vmin=0, vmax=25)

# Plot the boundary of the Natural Park
gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=2)

# Plot the stations
ax.plot(sta_lon, sta_lat, 'k^', mec='white')

# Customize the colorbar
cbar = plt.colorbar(em, ax=ax)

# Set custom ticks and labels for the colorbar
cbar.set_ticks([0, 12.5, 25])  # Ticks at Low (0), Medium, and High
cbar.set_ticklabels(['High (Error 0 %)', 'Medium (Error 12.5 %)', 'Low (Error 25 %)'])
# Move the colorbar label to the top
cbar.ax.set_title("Confidence")
# Labels and title
ax.set_xlabel("Longitude (º)")
ax.set_ylabel("Latitude (º)")

# Save the figure
plt.savefig('Error2.png', dpi=600)
plt.savefig('Error2.pdf', dpi=600, format='PDF')

# Show the plot
plt.show()
