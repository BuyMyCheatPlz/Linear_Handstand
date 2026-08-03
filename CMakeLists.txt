cmake_minimum_required(VERSION 3.14)
project(LQR_InvertedPendulum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Eigen3 ──────────────────────────────────────
# 方式 1: 系统安装的 Eigen (推荐)
#   Ubuntu:  sudo apt install libeigen3-dev
#   macOS:   brew install eigen
#   Windows: vcpkg install eigen3
find_package(Eigen3 QUIET)

if(NOT Eigen3_FOUND)
    # 方式 2: 手动指定 Eigen 路径
    #   修改下行路径为你的 Eigen 目录（包含 Eigen/ 子目录的父目录）
    #   例如: set(EIGEN3_INCLUDE_DIR "C:/libs/eigen-3.4.0")
    set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Eigen include directory")

    if(EIGEN3_INCLUDE_DIR)
        include_directories(${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR
            "Eigen3 未找到！\n"
            "  请安装 Eigen3，或在 CMake 中设置 EIGEN3_INCLUDE_DIR 指向 Eigen 的 include 目录。\n"
            "  例如: cmake -DEIGEN3_INCLUDE_DIR=/path/to/eigen ..\n"
            "  下载: https://eigen.tuxfamily.org/")
    endif()
else()
    include_directories(${EIGEN3_INCLUDE_DIRS})
endif()

# ── 可执行目标 ──────────────────────────────────
add_executable(lqr_pendulum
    main.cpp
    lqr_inverted_pendulum.cpp
    lqr_inverted_pendulum.h
)

# 平台特定设置
if(MSVC)
    # MSVC: 允许 M_PI 等数学常量
    target_compile_definitions(lqr_pendulum PRIVATE _USE_MATH_DEFINES)
endif()

# GCC 特定: 部分工具链 (如某些 MSYS2/MinGW) 在 -O1/-O2 下有
# std::ofstream ABI 问题，Debug 模式下默认 -O0 安全。
# Release 模式下可用 -Og 做轻度优化。
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(lqr_pendulum PRIVATE -Wall -Wextra)
endif()
