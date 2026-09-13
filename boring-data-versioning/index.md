---
title: LFS should be boring
published: 2026-09-13
author: David Muhr
abstract: "Large-file versioning should be simple and fast: no extra service, no separate workflow, no pipeline system. Just files that stay aligned with the commits that reference them."
---

I want a fairly unremarkable thing from data versioning: check out an old Git commit and get the files that belong to it. Not the latest dataset. Not the model someone happens to have on their laptop. The files that belong to that commit.

The usual answers are Git LFS or DVC. Both are reasonable, and both solve substantially more than a shell script that uploads a directory. But they make different decisions about the scope of data storage and where the additional machinery should live.

I think those decisions are more interesting than the usual feature comparison.^[Disclosure: I built [Gat](https://getgat.dev/), a large-file versioning tool.]

## Git LFS, the cost of another service

Git LFS fits closely into Git's existing workflow. Its clean filter replaces file contents with a pointer when adding them to Git; its smudge filter restores the contents when checking them out. A pre-push hook uploads the required LFS objects before publishing the commits.^[These operations are described in the [Git LFS specification](https://github.com/git-lfs/git-lfs/blob/main/docs/spec.md).]

That integration is useful. Once configured, the familiar Git commands do much of the work.
The trade-off is that **restoring a working tree now depends on an additional service**.

The normal LFS setup speaks to an LFS endpoint, not directly to an arbitrary S3 bucket. A hosting provider can operate that endpoint, or it can be self-hosted. Standalone custom transfer agents can also bypass the API server, so saying that LFS *always requires a dedicated server* would be inaccurate. But using storage you already control is not necessarily just a matter of supplying its bucket URL.^[Git LFS documents [custom transfer agents](https://github.com/git-lfs/git-lfs/blob/main/docs/custom-transfers.md), including operation without an API server.]

For a project whose Git host already provides suitable LFS storage, this can be an excellent arrangement. However, I would like to be able to choose whatever storage is practical without an additional service necessary to retrieve data.

## DVC and the cost of generality

DVC addresses that storage requirement directly. It can keep data in S3, Azure Blob, GCS, or a local directory, while Git versions the metadata. It also supports pipelines, parameters, metrics, and experiments.^[See DVC's [remote-storage documentation](https://doc.dvc.org/user-guide/data-management/remote-storage) and [user guide](https://doc.dvc.org/user-guide).]

There is an important qualification here: **none of those pipeline features is required for basic file versioning**. A workflow built around `dvc add`, `dvc push`, and `dvc pull` is supported. However, it appears that a good portion of the engineering effort has been put into the larger feature set instead of making git-based large file storage the main priority.

For me, the question is one of scope. I do not need a tool that can also describe how the data was produced. I only need to record which bytes belong to which commit. DVC can do that; but its apparant that it's not the only thing it wants to do.

## A deliberately boring alternative

That narrower responsibility is what I am trying to keep in [Gat](https://github.com/getgat-dev/gat): version large files with Git, leave the rest of the workflow alone. Essentially, I want a **lock file and a highly-efficient engine to process data around it**.

Git commits `gat.lock`, which maps paths to content IDs. The bytes live in a cache and, when configured, an object store or a directory. `gat add` records content, `gat push` uploads it, and `gat pull` downloads and restores it.

The entire engineering effort is deliberately spent on that problem. It is written in Rust, with a clear focus on keeping operations blazingly fast: massively parallel hashing, state-aware synchronization that avoids unnecessary rehashing, copy-on-write materialization to avoid copying cached files and many more optimizations.^[See the [Gat repository](https://github.com/getgat-dev/gat) and [performance guide](https://getgat.dev/guides/improving-performance). Reflinks require filesystem support and configuration; the default materialization strategy is an independent copy.]

I am not suggesting replacing an LFS or DVC setup that already works just because another tool is available. The difference I am aiming for is solving one problem well such that data versioning in Git becomes less of a chore.
