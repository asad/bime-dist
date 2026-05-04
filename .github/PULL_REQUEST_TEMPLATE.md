<!--
Thank you for your contribution to BIME!

Please:
  - Read CONTRIBUTING.md before opening this PR.
  - Ensure `node tools/run-tests.js` reports 0 failures locally.
  - Add a CHANGELOG entry under [Unreleased] describing the change.
  - Keep the PR small and focused — one logical change.
-->

## Summary

<!-- One-sentence description of what this PR does. -->

## Type of change

<!-- Mark with [x]. -->

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would change existing behaviour)
- [ ] Documentation update
- [ ] Test-only change
- [ ] Refactor / internal change (no API change)
- [ ] Performance improvement
- [ ] New molecule(s) added to `common-molecules.js`

## Linked issue(s)

<!-- e.g. "Closes #123" or "Refs #45" -->

## What changed and why

<!-- Bullet-point summary of the changes and the motivation. -->

## How was this tested?

- [ ] `node tools/run-tests.js` — all tests pass
- [ ] `node tools/audit-aromatic.js` — 0 issues across all molecules
- [ ] `node tools/build.js` — bundle rebuilds, bundle tests pass
- [ ] Manually verified in Chrome
- [ ] Manually verified in Firefox
- [ ] Manually verified in Safari
- [ ] Verified on touch device (iPad / Android)
- [ ] Other: <!-- describe -->

## Screenshots (UI changes only)

<!-- Drag-drop before/after images here. -->

## Checklist

- [ ] My code follows the project style (ES5-only in `editor/*.js`, no
      external runtime dependencies, strict equality).
- [ ] My code is CSP-safe (no `eval`, no `new Function`, no unsafe
      `innerHTML` on user data).
- [ ] I have added tests that prove my fix is effective or that my
      feature works.
- [ ] I have updated `docs/USAGE.md` if I changed user-visible behaviour.
- [ ] I have updated `CHANGELOG.md` under `[Unreleased]`.
- [ ] My commits are signed with my real name and a working email.
- [ ] I am the author of every line of code in this PR (or have the
      right to submit it under Apache 2.0).
- [ ] I have NOT introduced any third-party code without proper
      attribution and licence compatibility.

## Notes for reviewers

<!-- Anything specific you want reviewers to focus on? -->
