import matplotlib.pyplot as plt

# Create a figure and axes
fig, ax = plt.subplots()

# Plot your data
x = [1, 2, 3, 4, 5]
y=x
#y = [2, 4, 6, 8, 10]
ax.plot(x, y)

# Enable interactive mode
plt.ion()

# Show the plot
plt.show()

# Track mouse click locations
while True:
    # Wait for a mouse click event
    clicks = plt.ginput(n=1, timeout=10)

    # Check if a mouse click occurred
    if clicks:
        # Get the x and y coordinates of the click
        x_click, y_click = clicks[0]

        # Print the click location
        print(f"Clicked at ({x_click}, {y_click})")
