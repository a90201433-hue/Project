CC = g++
CFLAGS = -g -c -Wall -Ilib
LDFLAGS =
SRCDIR = src
OBJDIR = obj
EXECUTABLE = test

# Находим все .cpp файлы в src/
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# Преобразуем пути src/*.cpp → build/*.o
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

# Создание директории build при необходимости
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Правило для сборки .o файлов в build/
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) $< -o $@

pict: $(EXECUTABLE)
	./$(EXECUTABLE)
	python3 picture.py

run: $(EXECUTABLE)
	./$(EXECUTABLE)
	python3 animation.py
	python3 picture.py
	python3 Err.py

clean:
	rm -rf $(OBJDIR) $(EXECUTABLE)

