# A Divide-and-Conquer Approach for Global Orientation of Non-Watertight Scene-Level Point Clouds Using 0-1 Integer Optimization(DACPO)

This paper is published at ACM Transactions on Graphics (SIGGRAPH 2025). Please refer to our [paper](https://arxiv.org/abs/2505.23469).

![1750858552307](image/README/1750858552307.png)

### Tested Platform

- Windows 11
- Visual Studio 2022
- 11th Gen Intel(R)Core(TM)i9-11900K 3.50GHz
- 64GB RAM

### Dependencies
- pcl
- open3d
- Eigen
- boost

**1. PCL**

Install using vcpkg.

`vcpkg install pcl:x64-windows`

`vcpkg integrate project`

After executing `vcpkg integrate project`, you will see a prompt:

> With a project open, go to Tools->NuGet Package Manager->Package Manager Console and paste:
> Install-Package "vcpkg.D.sdk.vcpkg" -Source "C:\Users\19892"

Follow the instructions as prompted.

**2. Open3D**

Download the Open3D package from [Github releases](https://github.com/isl-org/Open3D/releases), and configure the header files and library files in the solution manager. Link open3d.lib in the linker.

The header files include:

```
include/
include/3dparty 
```
Place the dynamic libraries in the executable file directory.



**3. Eigen**

Install using vcpkg:

`vcpkg install Eigen3:x64-windows`

`vcpkg integrate project`



**4. Boost**

Install using vcpkg:

`vcpkg install boost:x64-windows`

`vcpkg integrate project`


**5. Configure Header and Library Paths**

In all configurations, set the VC++ directory include paths:

```
include:
src/include
src/include/3dparty
```

### Build with CMake

If you prefer to use CMake instead of Visual Studio, you can use the provided CMakeLists.txt:

1. **Create build directory:**
```bash
mkdir build
cd build
```

2. **Configure the project:**
```bash
# Using vcpkg toolchain (recommended)
cmake .. -DCMAKE_TOOLCHAIN_FILE=[path_to_vcpkg]/scripts/buildsystems/vcpkg.cmake

# Or use system default
cmake ..
```


### Dataset
Our code can process a single scene in a single ply file.
Batch preprocess dataset and batch running code will be released later.

#### ScanNet
please refer to [ScanNet](https://github.com/ScanNet/ScanNet)

#### SceneNN
please refer to [SceneNN](https://github.com/zhaoyu-zhao/SceneNN)




### Citation

If you find this work useful in your research, please consider citing:

```bibtex
@misc{li2025divideandconquerapproachglobalorientation,
      title={A Divide-and-Conquer Approach for Global Orientation of Non-Watertight Scene-Level Point Clouds Using 0-1 Integer Optimization}, 
      author={Zhuodong Li and Fei Hou and Wencheng Wang and Xuequan Lu and Ying He},
      year={2025},
      eprint={2505.23469},
      archivePrefix={arXiv},
      primaryClass={cs.CV},
      url={https://arxiv.org/abs/2505.23469}, 
}
```

