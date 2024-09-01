import matplotlib.pyplot as plt

# Список файлов для обработки
file_list = ["analytical_for_cubic.txt", "cubic.txt", "analytical_for_linear.txt", "linear.txt"]

# Создаем общее окно для всех графиков с 2 столбцами и количеством строк, необходимым для всех графиков
num_files = len(file_list)
num_columns = 2
num_rows = (num_files + 1) // num_columns

fig, axes = plt.subplots(num_rows, num_columns, figsize=(12, 8))

# Для каждого файла в списке
for i, filename in enumerate(file_list):
    # Инициализация списков для хранения данных
    t_values = []
    u_values = []

    # Открываем файл для чтения
    with open(filename, 'r') as file:
        # Считываем данные из файла
        for line in file:
            # Разделяем строку на значения t и u
            t, u = map(float, line.strip().split())
            t_values.append(t)
            u_values.append(u)

    # Вычисляем текущую позицию в сетке подграфиков
    row = i // num_columns
    col = i % num_columns

    # Строим график в соответствующем подграфике
    ax = axes[row, col]
    ax.plot(t_values, u_values)
    ax.set_xlabel('t')
    ax.set_ylabel('u')
    ax.set_title(f'Plot for {filename}')

# Устанавливаем расстояние между подграфиками и скрываем пустые подграфики
for i in range(num_files, num_rows * num_columns):
    row = i // num_columns
    col = i % num_columns
    fig.delaxes(axes[row, col])

# Отображаем все графики в одном окне
plt.tight_layout()
plt.show()
