import matplotlib.pyplot as plt

# Открываем файл для чтения
with open("data.txt", "r") as file:
    lines = file.readlines()

# Создаем списки для хранения данных
x = []
y_values = []

# Считываем данные из файла
for line in lines:
    data = line.strip().split()
    x.append(float(data[0]))
    y_values.append([float(val) for val in data[1:]])

# Определение количества y-значений
num_y_values = len(y_values[0])

# Создаем графики для каждого y с отдельными осями y
fig, axs = plt.subplots(num_y_values, 1, sharex=False, figsize=(8, 6))

# Настраиваем каждый график
for i, ax in enumerate(axs):
    
    y_data = [y[i] for y in y_values]
    ax.plot(x, y_data, label=f"x{i+1}")
    ax.set_ylabel(f"y{i+1}")

# Настраиваем оси x для последнего графика
axs[-1].set_xlabel("t")

# Отображаем графики
plt.show()
