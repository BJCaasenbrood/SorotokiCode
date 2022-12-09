# Installation

The toolkit is easy to install. The fastest and easiest way to acquire the toolkit is through [git](https://git-scm.com/downloads). We highly recommend using Git, since this allows you also to update the toolkit. 

???+ warning "Prerequisite toolboxes of MATLAB"
    To use Sorotoki, make sure you have the following dependencies installed:

    - [Optimization Toolbox](https://nl.mathworks.com/products/optimization.html) 
    - [Partial Differential Equation Toolbox](https://nl.mathworks.com/products/pde.html)
    - [Image Processing Toolbox](https://nl.mathworks.com/products/image.html) 
    - [Matlab Coder](https://www.mathworks.com/products/matlab-coder.html)


## Step-by-step installation  
##### Step 1: Downloading the Sorotoki installation folder
You can directly clone the repository using the command:

``` git
git clone --depth 1  https://github.com/BJCaasenbrood/SorotokiCode.git
```

Alternatively, you can download the latest version of the toolkit below, and unpack the compressed folder at any desired working directory. We do recommend that the installation folder of the toolkit named: `SorotokiCode` is under the folder `...\Documents\Matlab\` (at least for Windows and Linux OS's).

- Toolkit installation size (only code):     `~20 MB` :feather:
- Toolkit installation size (with STL data): `~120 MB` 🗿

??? Direct-Download
    [Sorotoki v.2.19.stable (code)](#){ .md-button .md-button--primary }
    [Sorotoki v.2.19.stable (code + STLs)](#){ .md-button .md-button--primary }


##### Step 3: Run the installation

Once all the prerequisites are (properly) installed, you first have to configure the toolkit with MATLAB's search paths. Setting up these paths is relatively straightforward. Simply run:

``` matlab
sorotoki();
```

The file should be located under `...\Documents\MATLAB\SorotokiCode\.` The toolkit will proceed the installation by building search paths in Matlab. The toolkit will request the installation of the prerequisite toolboxes if they are not installed. These include: Optimization Toolbox, Partial Differential Equation Toolbox, Image Processing Toolbox and Matlab Coder. Especially the latter ensures that embedded Sorotoki functions can be converted to `c`or `c++` equivalent code in the form of mex files; that greatly enhances computational performance. Finally, a verification routine is performed that checks if the toolkit is installed correctly.

##### Step 4: That's it folks
And that's it, Sorotoki is now ready to use. The `sorotoki` installation command can also be used to run demo, check updates, see version number, and find the documentation file. You can do these by running the following commands:

``` matlab
sorotoki demo       % list demos
sorotoki d          % ...

sorotoki version    % version (and do i need to update?)
sorotoki v          % ...

sorotoki doc        % documentation link
sorotoki i          % ...

```

## Darn', problems during installation?
Please let me know by:

1. Post an issue at: [https://github.com/BJCaasenbrood/SorotokiCode/issues](https://github.com/BJCaasenbrood/SorotokiCode/issues)