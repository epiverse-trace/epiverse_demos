# epiverse_demos

A project to demonstrate how to solve outbreak analytics tasks
using open-source packages from the 
[Epiverse-TRACE](https://epiverse-trace.github.io/),
[Epiforecasts](https://epiforecasts.io/), 
[Reconverse](https://www.reconverse.org/), and
[R Epidemics Consortium](https://www.repidemicsconsortium.org/projects/)
toolkits.

> [!NOTE]
> These entries are now hosted and maintained in the
> **How-to Guides** repository of the Epiverse-TRACE organisation.
> You can access them at: <https://epiverse-trace.github.io/howto/>

## Usage

1. Open any `task-...` R file in root.
2. Run the content to explore how a set of task are coded to solve an end goal.

All data is simulated or stored within R packages.


## Posit Cloud

The Posit Cloud already have all packages and dependencies installed. 

To keep record of your edits or notes,
Click on **Save a Permanent Copy**.


## Support

Join the [Discussions on GitHub](https://github.com/orgs/epiverse-trace/discussions).


## License

[MIT](https://choosealicense.com/licenses/mit/)


## Contributing

Contributions are always welcome!

In the [GitHub repository](https://github.com/epiverse-trace/epiverse_demos), 
fill issues or fork the repository to create a Pull Request. 


## Installation

1. Fork and clone the repository: <https://github.com/epiverse-trace/epiverse_demos>
2. Restore the R environment running:

```r
if(!require("renv")) install.packages("renv")
renv::restore()
```

