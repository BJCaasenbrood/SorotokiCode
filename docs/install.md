---
layout: default
title: Installation 
nav_order: 2
---

# Installation
{: .no_toc }

The fastest and easiest way to acquire the toolkit is through [**git**](https://git-scm.com/downloads). You can directly clone the repository using the git command:

```fortran
 git clone https://github.com/BJCaasenbrood/SorotokiCode.git
```

Alternatively, you can download the latest stable version below, and unpack the compressed folder at any desired directory. 

[Stable V3.03 (.zip)](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master){: .btn .btn-green .fs-5 .mb-4 .mb-md-0 .mr-2} [Stable V3.03 (.tar)](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master){: .btn .btn-green .fs-5 .mb-4 .mb-md-0 .mr-2} [View on Github](https://github.com/BJCaasenbrood/SorotokiCode){: .btn .fs-5 .mb-4 .mb-md-0}  


## Set-up and getting started
Before we get started, we first have to configure the toolkit with MATLAB's function paths. Setting up the toolkit is relatively straightforward. To set-up the toolkit, simply run the command below. That's it, the soft robotics toolkit is now ready-to-use. 
```rust
% installation command
sorotoki();
```
The command above is also used to update the toolkit. It is recommended to run `sorotoki.m`{: .text-purple-000} to check for updates occasionally. 

To get started, type the following line in the command window:
```rust
% show demos
sorotoki('demo');
```

