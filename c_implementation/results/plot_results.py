import matplotlib.pyplot as plt


def plot_points(x_points, y_points, labels):
    unique_labels = list(set(labels))

    label_to_color = {label: i for i, label in enumerate(unique_labels)}

    colors = [label_to_color[label] for label in labels]

    plt.scatter(x_points, y_points, c=colors, cmap='viridis')  # Scatter plot with colors based on labels
    plt.xlabel('x')  # x-axis label
    plt.ylabel('y')  # y-axis label
    plt.title('Plot of x and y points with colors based on labels')  # Title of the plot
    plt.grid(True)  # Show grid
    # plt.colorbar(label='Labels')  # Add colorbar with label names
    plt.show()  # Display the plot


def read_points_from_file(path):
    x_points = []
    y_points = []
    labels = []

    with open(path, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            x_points.append(float(values[0]))
            y_points.append(float(values[1]))
            labels.append(int(values[2]))

    return x_points, y_points, labels


if __name__ == '__main__':
    file_path = 'results_sequential.txt'
    x, y, lab = read_points_from_file(file_path)
    plot_points(x, y, lab)
