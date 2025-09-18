cmake_minimum_required(VERSION 3.10)
project(TwoFluidSimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configurar directorios de inclusión
include_directories(include)

# Archivos fuente
set(SOURCES
    src/main.cpp
    src/dispersion.cpp
    src/simulator.cpp
)

add_executable(TwoFluidSimulator ${SOURCES})

# Habilitar optimizaciones
target_compile_options(TwoFluidSimulator PRIVATE -O3 -march=native)

# Buscar todos los .py en la carpeta raíz del proyecto
file(GLOB PY_SCRIPTS "${CMAKE_SOURCE_DIR}/*.py")

# Copiar cada uno después del build
foreach(script ${PY_SCRIPTS})
    add_custom_command(TARGET TwoFluidSimulator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${script}"
        $<TARGET_FILE_DIR:TwoFluidSimulator>)
endforeach()

add_custom_command(TARGET TwoFluidSimulator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:TwoFluidSimulator>/data
)
