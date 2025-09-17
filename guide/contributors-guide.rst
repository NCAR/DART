.. index:: contributing, contribute

.. _contributors-guide:

Contributors Guide
===================

Welcome to the DART contributors guide! We appreciate your interest in 
contributing to the project. As an open-source project, we rely on the 
creativity and expertise of people like you. Whether you're fixing bugs, 
adding new features, improving documentation, or writing tests, please 
follow the contributors guide.

What Can I Do?
----------------

Contributors come from many backgrounds: scientists working with models, 
students learning data assimilation, software developers, general 
researchers and others. No matter your expertise, there are ways you can 
contribute.

* **Report Bugs**: If you find a bug, please report it by opening an issue on our 
  `GitHub repository <https://github.com/NCAR/DART/issues>`__.

* **Fix Bugs**: Look through the `issue tracker <https://github.com/NCAR/DART/issues>`__
  for bugs that need fixing. Please comment on the issue to let us know you're working on it
  or reach out to the DAReS team at dart@ucar.edu. This helps make sure your work isn't
  duplicating someone else's efforts.

* **Add Features**: If you have an idea for a new feature, open an issue to 
  discuss it before starting work.

* **Improve Documentation**: Help us improve our documentation by making it
  clearer and more comprehensive, or by simplifying complex concepts.

* **Write Tests**: Help us ensure the quality of DART by writing tests or 
  supplying test cases for your model.

Source Code
------------
The source code for DART is available on our GitHub repository `NCAR/DART <https://github.com/NCAR/DART>`__.
Feel free to explore and understand the codebase.


Reporting a Bug
----------------

If you find a bug, please report it by opening an issue on our `GitHub repository <https://github.com/NCAR/DART/issues>`__.
Include as much detail as necessary to help us understand and reproduce the issue.

A bug report should contain the following:

* The steps someone needs to take to reproduce the bug.
* What you expected to happen.
* What actually happened.
* Details about your system environment (e.g., compiler version, operating system, 
  MPI/NetCDF library version, etc). 
* Any relevant logs, error messages, or (if possible) a minimal reproducible example. 


Pull Requests
--------------

We welcome pull requests! Please take a read through this contributors guide
before developing with DART. As a small team maintaining a community code, 
we may not always be able to accept every contribution, but we greatly 
appreciate your efforts and interest.

DART follows the `GitHub Flow <https://guides.github.com/introduction/flow/>`__ workflow.

* **Fork the Repository**: Click the “Fork” button on the `NCAR/DART GitHub page <https://github.com/NCAR/DART>`__.

* **Clone Your Fork**: Clone your forked repository to your local machine.

   .. code-block:: bash

       git clone https://github.com/your-github-username/DART.git

* **Create a Branch**: Create a new branch for your work.

   .. code-block:: bash

       cd DART
       git checkout -b my-feature-branch

  Replace my-feature-branch with a name for your branch that describes your feature or fix.

* **Make Changes**: Make your changes and commit them.

  Be as helpful to your reviewers and the DAReS team as you can. Only include changes relevant
  to your issue and avoid changes that are not relevant to the particular issue. Keep
  PRs small and focused as this makes them much easier to review. 

* **Push to Your Fork**: Push your changes to your forked repository.

   .. code-block:: bash

       git push -u origin my-feature-branch

* **Open a Pull Request**: From your fork, open a pull request to the main repository. 
  
  Include a short description of your changes —the “why”— and reference any related issues.

* **Be responsive to the review**: Your pull request will be reviewed by the DAReS team.

  We may ask for changes or clarifications.
  Please be patient as we work through the review process. If accepted, your pull request will be
  merged into the main repository and released in the next version of DART.
