Contributors' guide
===================

Contributing to DART
--------------------

This section describes how you can contribute your work to DART. Because DART
is an open-source project, your contributions are welcome. Many user-provided
contributions have widely benefited the earth science community.

To ensure you aren't duplicating efforts, contact DAReS staff by emailing 
dart@ucar.edu before you expend considerable development time.

All of the source code is hosted in the `DART GitHub repository
<https://github.com/NCAR/DART>`__.

Before you start developing, you should be familiar with the `GitHub
workflow <https://guides.github.com/introduction/flow/>`_. The GitHub worflow 
involves:

1. Creating a *fork* of the DART project. A fork is a publically visible copy
   of the repository that is stored in your GitHub account.
2. Creating a *branch* for your feature with an appropriate name for your
   project, and when you are finished with your changes you can commit them
   to your fork. After testing locally on your machine, you can push them to
   your fork.

   .. Important::
   
      At this point, everyone can see the changes you made on your fork.

   When you are ready to begin the conversation about merging your work into
   the original project (called the DART repository master), you can create a
   pull request, which will show your changes. After reviewing and testing
   your changes, the pull request will be addressed appropriately by the DART
   development team.

Keeping your work private until you publish
-------------------------------------------

You may want to keep your work private until it is ready for publication or
public viewing.

Follow these steps to hide sensitive code until you are ready to contribute it 
to DART your work has been published.

1. First, create a public fork of the DART repository by following the 
   steps listed above.
2. Next, `create a private <https://help.github.com/en/articles/create-a-repo>`__
   repository on GitHub.com. The name of your private repository is arbitrary,
   since only you and your private collaborators can see it.
3. Add your public fork as a
   `remote repository <https://help.github.com/en/articles/adding-a-remote>`__
   of your private repository. Your remote repository can be named
   "public_fork" or "upstream."
4. Add additional team members, if necessary.
5. Instead of pulling and pushing from your public fork, `develop on your
   private repository <https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes>`__.

.. note::
   
   Only three collaborators are allowed on a free non-institutional private
   repository. DAReS staff can collaborate with you on your private repository,
   but keep this three collaborator limit in mind if you using a free GitHub
   account.
