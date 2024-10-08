---
title: The importance of local development
published: 2024-10-06
author: David Muhr
abstract: Local development is crucial for efficient software development, it enables faster iteration, better debugging and ensures consistency between local and production environments. 
--- 

A good developer experience (DX) increases productivity, reduces cognitive load, and improves developer satisfaction, 
leading to higher code quality, faster onboarding, and talent retention.
An easy-to-use local development environment might be the foundation for good DX.

I regularly think back to the
[hilarious discussion, where Elon Musk proposed to rewrite Twitter](https://www.youtube.com/watch?v=FkNkSQ42jg4) reminding
me about the following statement.

> One of the biggest problems is that you can't run Twitter locally. When I was at Facebook, you could run all of Facebook on a laptop. There is no way to run Twitter outside of a very bespoke configuration in a datacenter and it makes it nearly impossible to build anything new on it. <cite>George Hotz</cite>

The importance of a good developer experience seems to correlate with the number of people working on a project.
As a lone fighter, it's quite simple to spin-up a bespoke setup, because its configuration is known to a single person.

I was a lone fighter myself for much of the research journey during my PhD, which might be the reason why I
previously neglected developer experience, except for my own experience.

It's common to think about technical debt when adding new tools or features, but that is possibly not the
entire story. The concept of ==DX-debt== is equally relevant, yet it is not sufficiently discussed, in my opinion.

Nowadays, there are options such as [Draft](https://github.com/Azure/draft/),
[Skaffold](https://github.com/GoogleContainerTools/skaffold), [Tilt](https://github.com/tilt-dev/tilt) or
[Garden](https://github.com/garden-io/garden) to spin-up entire clusters locally, or provide
[pseudo-local](https://github.com/readme/guides/developer-onboarding) development environments.
It's never been easier to provide a great developer experience. What do you think?
