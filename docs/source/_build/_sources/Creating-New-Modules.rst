################################################
Creating-New-Modules
################################################


#. 
   **Create a New Module Directory**\ :


   * Navigate to the ``modules`` directory in the Kleborate project.
   * Create a new directory for your module. The directory name should be descriptive of the module's functionality.

#. 
   **Create Module Python File**\ :


   * Inside the newly created directory, create a Python file with the same name as the directory. This file will contain the module's implementation.
   * For example, if your module is named ``enterobacterales__species``\ , create a file named ``enterobacterales__species.py`` inside the ``enterobacterales__species`` directory.

#. 
   **Implement Module Functionality**\ :


   * Define the functionality of your module inside the Python file.
   * You need to define functions for adding CLI options, checking CLI options, checking external program dependencies, getting module headers, getting module results, and any other functionality your module requires.
   * Ensure that your module follows the conventions and guidelines of existing Kleborate modules.

#. 
   **Add CLI Options**\ :


   * Implement a function to add command-line options specific to your module. This function should accept an ``argparse.ArgumentParser`` object as an argument and add options using its methods.
   * The function should be named ``add_cli_options(parser)`` and return the argument group created for the module's options.

#. 
   **Check CLI Options**\ :


   * Implement a function to check the CLI options provided by the user. Ensure that the provided options are valid and meet the module's requirements.
   * This function should be named ``check_cli_options(args)`` and raise appropriate errors or warnings if the options are invalid.

#. 
   **Check External Program Dependencies**\ :


   * Implement a function to check the external program dependencies required by your module.
   * This function should be named ``check_external_programs()`` and return a set containing the names of external programs required by the module.

#. 
   **Define Prerequisite Modules**\ :


   * Define a function named ``prerequisite_modules()`` that returns a list of module names that your module depends on.
   * This ensures that prerequisite modules are executed before your module.

#. 
   **Define Module Headers**\ :


   * Implement a function to define the headers for the module's output.
   * This function should return two lists: one for the full headers and one for the stdout headers.

#. 
   **Define Module Results**\ :


   * Implement a function to get the results produced by your module.
   * This function should accept necessary arguments like assembly, minimap2 index, command-line arguments, and other required data.
   * It should return a dictionary containing the results.

#. 
   **Test Your Module**\ :


   * Test your module with different inputs and scenarios to ensure that it behaves as expected.
   * Verify that the module produces correct results and handles errors.

#. 
   **Document Your Module**\ :


   * Write documentation for your module, including its purpose, usage, options, and any other relevant information.
   * Update the main Kleborate documentation to include information about your module.
