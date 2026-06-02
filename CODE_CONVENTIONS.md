# Code conventions

openLowdin has been developed since 2008 by many different young scientific developers. We have been trying to follow some degree of standardisation in the code, but we are still far from where we would like to be due to pressure from our academic duties and, let's be honest, lack of discipline. 

In any case, in this document, we are listing the desired openlowdin code conventions as a guideline for any developer. Please follow these instructions for your new implementations and help us to bring order if you find that any of these conventions are not followed in the code.  

## Code Structure

openLowdin depends on three types of code:

- in-house code developed by the team and stored in `src/` folder 
- external code developed by third parties but distributed within openlowdin in the `utilities` folder for simplicity
    - open-source code (we don't want legal issues)
    - not so common software, not available in package managers (e.g. apt, pip, dnf)
    - slightly adapted code to work with openlowdin
- external software dependencies not distributed with openlowdin
    - commonly available in package managers (e.g. apt, pip, dnf) 
    - large code, long compilation times

## Code development workflow 

For major implementations:

1. Please contact the team! <https://github.com/efposadac/openLOWDIN> or <https://openlowdin.github.io/openLOWDIN_manual/>, we can guide you and discuss if the proposed changes can be included in the main repository, in some cases, a scientific collaboration could be established. 
2. Check what has been implemented already in `src` let's avoid code repetitions. 
3. Create your own git branch. 
4. Try to write modular code, avoid spaghetti code between libraries, and follow the global structure:
    ```
    src/lowdin.x.f90  ( main program, read inputs, construct objects, call Solver )
    src/Solver/  ( library, call the requested methods )
    src/core/  ( library, contains general modules for loading input, building main objects, general math functions )
    src/libraryname/  ( library, specific for given task/method, e.g. computing integrals, HF, DFT )
    ```
5. Follow Fortran standard section above
6. Add any external dependencies into the `configure` file or `utilities` folder
7. The basics, compile and link: `make` and `make install` 
8. Create your own test! make sure that the new features are used in a calculation.  
    - Create a `.lowdin` input in `test` folder
    - Create a `.py` script in `test` folder to automatically compare with your reference value. You can follow any of the tests there as an example, notice the common functions to extract information from the output are include in `test/lowdinTestFunctions.py` remember, avoid code repetitions
9. Run `make test` all tests should pass! not failed tests, no code crashs, no Balrogs!
10. `git commit`, `git push`, `git merge` and let's wait for the team to approve merge. 
11. Document your changes in the manual <https://openlowdin.github.io/openLOWDIN_manual/>


## Fortran Standard

1. Follow Fortran 2003 standard, other programming languages are welcome but... we like Fortran. 
2. Use English for comments and variable naming. Yes, we know, all the original developers are Spanish speakers, but English is the scientific language. Except for funny comments! 
3. Module name convention: First capitalized letter, name of the library (optional), lower case unique name. 
4. Subroutine/function name convention: Name of the module, `"_"` symbol, lower case unique name, use upper case letter to mark the beginning of a new word (e.g. `Mymodule_myNewFunction`) 
5. Public type name convention: Name of the module, `"_"` symbol, lower case unique name.
6. General name convention: Lower case unique name, use upper case letter to mark the beginning of a new word (e.g. `myNewVariable`), use `"_"` symbol if necessary for readability or for following a common convention
7. Formatting: Indentation of sections of code should use two spaces, not tabs. We used `fprettify --indent 2` for auto-formatting Fortran files <https://github.com/fortran-lang/fprettify>.
8. Add comments to explain code in simple words. 
9. We have a mess with data types... but 64 bits should be the default for everything, except for optimization cases in really large arrays
