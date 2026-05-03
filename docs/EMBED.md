# Embedding BIME

This guide explains how to drop BIME into your own webpage, blog post,
LMS course, dashboard, or interactive notebook.

BIME is a pure-JavaScript library — no React, no Vue, no Angular, no
build step. You can embed it in three lines of HTML.

---

## Table of contents

1. [The minimal embed (3 lines of HTML)](#1-the-minimal-embed-3-lines-of-html)
2. [Iframe embed (zero-touch)](#2-iframe-embed-zero-touch)
3. [Production embed with the bundle + SRI](#3-production-embed-with-the-bundle--sri)
4. [Embedding inside an LMS](#4-embedding-inside-an-lms)
5. [Embedding in a Jupyter / JupyterLab notebook](#5-embedding-in-a-jupyter--jupyterlab-notebook)
6. [Embedding in Observable](#6-embedding-in-observable)
7. [Embedding inside a Markdown blog (Hugo / Jekyll / 11ty / Astro)](#7-embedding-inside-a-markdown-blog-hugo--jekyll--11ty--astro)
8. [Embedding in React, Vue, Svelte (without rebuilding BIME)](#8-embedding-in-react-vue-svelte-without-rebuilding-bime)
9. [The JavaScript API for embedded use](#9-the-javascript-api-for-embedded-use)
10. [Capturing student work from a parent page](#10-capturing-student-work-from-a-parent-page)
11. [CSP and security considerations](#11-csp-and-security-considerations)

---

## 1. The minimal embed (3 lines of HTML)

This is all you need:

```html
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="https://asad.github.io/bime-dist/css/style.css">
</head>
<body>

    <!-- 1. Container -->
    <div id="my-editor" style="width:100%;height:500px"></div>

    <!-- 2. BIME source modules -->
    <script src="https://asad.github.io/bime-dist/editor/Molecule.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/Layout.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/Templates.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/SmilesParser.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/SmilesWriter.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/MolfileWriter.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/Renderer.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/History.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/Tools.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/CipStereo.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/SmartsParser.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/SmartsMatch.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/SmartsWriter.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/ImageExport.js"></script>
    <script src="https://asad.github.io/bime-dist/editor/MolEditor.js"></script>

    <!-- 3. Initialize -->
    <script>
        var editor = new MolEditor('my-editor', '100%', '500px');
        editor.readGenericMolecularInput('c1ccccc1');   // load benzene
    </script>

</body>
</html>
```

Open the file in a browser — you have an embedded BIME editor.

---

## 2. Iframe embed (zero-touch)

If you can't add `<script>` tags (e.g. in some restrictive blogging
platforms), use an iframe pointing to the live workbench:

```html
<iframe src="https://asad.github.io/bime-dist/workbench.html"
        width="100%"
        height="700"
        title="BIME molecular editor"
        loading="lazy"
        referrerpolicy="no-referrer">
</iframe>
```

This embeds the **full workbench** (toolbar, sidebar, examples). If you
want a stripped-down editor-only view, use the script-tag embed in §1.

---

## 3. Production embed with the bundle + SRI

For a smaller HTTP footprint and tamper-detection, use the pre-built
single-file bundle from `dist/`:

```html
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="https://asad.github.io/bime-dist/css/style.css">
</head>
<body>

    <div id="editor" style="width:100%;height:500px"></div>

    <!-- Single bundled file with Subresource Integrity -->
    <script src="https://asad.github.io/bime-dist/dist/bime.min.js"
            integrity="sha384-<paste-from-dist/SRI.txt>"
            crossorigin="anonymous"></script>

    <script>
        var editor = new MolEditor('editor', '100%', '500px');
        editor.readGenericMolecularInput('CCO');
    </script>

</body>
</html>
```

Get the current `integrity` value from
[`dist/SRI.txt`](https://github.com/asad/bime-dist/blob/main/dist/SRI.txt) —
it changes with every release.

The bundle is ~480 kB minified, ~120 kB gzipped — small enough to load
without delay over typical broadband.

---

## 4. Embedding inside an LMS

### Moodle

1. **Add an HTML block** to your course page.
2. **Switch the editor to source view** (`<>` icon).
3. **Paste an iframe** from §2 above.
4. Save.

Your students see the workbench inside the course page, with no login
or navigation away.

### Canvas

1. **Edit a page**.
2. **Switch to HTML view**.
3. Paste the iframe.
4. Save.

Canvas may sandbox the iframe; if it does, make sure the URL is added
to your institution's allowlist.

### Blackboard

1. **Edit a content area → Build Content → HTML**.
2. Paste the iframe.
3. Save.

### Google Classroom

Google Classroom doesn't support iframe embedding, but you can post a
**link** to the BIME workbench (live demo or your self-hosted copy)
as an assignment. Students click and the workbench opens in a new tab.

For a guided tutorial flow within the LMS, see the next section on
[capturing student work](#10-capturing-student-work-from-a-parent-page).

---

## 5. Embedding in a Jupyter / JupyterLab notebook

Use a `%%html` cell:

```python
%%html
<div id="bime-cell" style="width:100%;height:500px"></div>
<script src="https://asad.github.io/bime-dist/dist/bime.min.js"></script>
<script>
  var ed = new MolEditor('bime-cell', '100%', '500px');
  ed.readGenericMolecularInput('Cn1c(=O)c2c(ncn2C)n(C)c1=O');  // caffeine
</script>
```

Or use IPython's `display(HTML(...))`:

```python
from IPython.display import HTML, display

display(HTML("""
<div id='bime-py' style='width:100%;height:500px'></div>
<script src='https://asad.github.io/bime-dist/dist/bime.min.js'></script>
<script>
  var ed = new MolEditor('bime-py', '100%', '500px');
  ed.readGenericMolecularInput('CCO');
</script>
"""))
```

For JupyterLite (browser-only Jupyter), the same works — BIME runs
client-side just like the rest of JupyterLite.

### Bidirectional Python ↔ BIME

To pass SMILES from Python to BIME and back, use `IPython.display` to
inject a button:

```python
%%html
<div id="bime-bd" style="width:100%;height:400px"></div>
<button id="grab">Grab SMILES</button>
<pre id="output"></pre>
<script src="https://asad.github.io/bime-dist/dist/bime.min.js"></script>
<script>
  var ed = new MolEditor('bime-bd', '100%', '400px');
  ed.readGenericMolecularInput('c1ccccc1');
  document.getElementById('grab').onclick = function() {
      document.getElementById('output').textContent = ed.getSMILES();
  };
</script>
```

---

## 6. Embedding in Observable

```js
{
    const div = html`<div style="width:100%;height:500px"></div>`;
    await require('https://asad.github.io/bime-dist/dist/bime.min.js');
    const editor = new MolEditor(div, '100%', '500px');
    editor.readGenericMolecularInput('Cn1c(=O)c2c(ncn2C)n(C)c1=O');
    return div;
}
```

Or use `htl` for strict CSP compliance — the BIME bundle is
require-compatible.

---

## 7. Embedding inside a Markdown blog (Hugo / Jekyll / 11ty / Astro)

Most static-site generators allow raw HTML in Markdown. Just paste:

```html
<div id="post-bime" style="width:100%;height:500px"></div>
<script src="https://asad.github.io/bime-dist/dist/bime.min.js"></script>
<script>
  new MolEditor('post-bime', '100%', '500px')
      .readGenericMolecularInput('OC(=O)CC(O)(CC(=O)O)C(=O)O');  // citric acid
</script>
```

If your generator escapes raw HTML, wrap it in a shortcode (Hugo) or a
`{% raw %}` block (Jekyll), or use an iframe (§2).

---

## 8. Embedding in React, Vue, Svelte (without rebuilding BIME)

You don't need to npm-install or rebuild BIME — load it once via a
script tag and call `MolEditor` from your framework's lifecycle:

### React

```jsx
import { useEffect, useRef } from 'react';

function BimeEditor({ initialSmiles = 'c1ccccc1', height = 500 }) {
    const ref = useRef(null);
    useEffect(() => {
        // Load BIME once if not already present
        if (!window.MolEditor) {
            const script = document.createElement('script');
            script.src = 'https://asad.github.io/bime-dist/dist/bime.min.js';
            script.onload = () => mount();
            document.body.appendChild(script);
        } else {
            mount();
        }
        function mount() {
            const editor = new window.MolEditor(ref.current, '100%', height + 'px');
            editor.readGenericMolecularInput(initialSmiles);
        }
    }, [initialSmiles, height]);
    return <div ref={ref} style={{ width: '100%', height }} />;
}

export default BimeEditor;
```

### Vue 3

```vue
<template>
  <div ref="container" :style="{ width: '100%', height: height + 'px' }"></div>
</template>

<script setup>
import { onMounted, ref } from 'vue';
const props = defineProps({
    initialSmiles: { type: String, default: 'c1ccccc1' },
    height:        { type: Number, default: 500 }
});
const container = ref(null);

onMounted(async () => {
    if (!window.MolEditor) {
        await new Promise((res, rej) => {
            const s = document.createElement('script');
            s.src = 'https://asad.github.io/bime-dist/dist/bime.min.js';
            s.onload = res; s.onerror = rej;
            document.body.appendChild(s);
        });
    }
    const editor = new window.MolEditor(container.value, '100%', props.height + 'px');
    editor.readGenericMolecularInput(props.initialSmiles);
});
</script>
```

### Svelte

```svelte
<script>
  import { onMount } from 'svelte';
  export let initialSmiles = 'c1ccccc1';
  export let height = 500;
  let container;

  onMount(async () => {
    if (!window.MolEditor) {
      await new Promise((res, rej) => {
        const s = document.createElement('script');
        s.src = 'https://asad.github.io/bime-dist/dist/bime.min.js';
        s.onload = res; s.onerror = rej;
        document.body.appendChild(s);
      });
    }
    const editor = new window.MolEditor(container, '100%', height + 'px');
    editor.readGenericMolecularInput(initialSmiles);
  });
</script>

<div bind:this={container} style="width:100%;height:{height}px"></div>
```

---

## 9. The JavaScript API for embedded use

Once you have a `MolEditor` instance:

```js
var editor = new MolEditor('container-id', '100%', '500px');

// Load
editor.readGenericMolecularInput('CCO');         // SMILES, MOL, or anything BIME recognises
editor.readSmiles('c1ccccc1');
editor.readMolfile(molV2000Text);

// Read
editor.getSMILES();          // canonical SMILES
editor.getCanonicalSmiles(); // alias
editor.getMolfile();          // MOL V2000
editor.getMolfileV3000();     // MOL V3000
editor.getMolecule();         // raw Molecule object — use SmilesWriter.write etc.

// Export
editor.exportSVG();           // SVG string
editor.exportPNG(2);          // PNG data URL at 2× resolution
editor.copySMILESToClipboard();
editor.copySVGToClipboard();

// Control
editor.clear();               // wipe canvas
editor.undo();
editor.redo();

// Callbacks for parent-page integration
editor.onAfterStructureModified = function(canonicalSmiles) {
    console.log('User changed structure → ' + canonicalSmiles);
};
editor.onAtomClicked = function(atom) {
    console.log('Atom clicked:', atom.symbol, atom.id);
};
editor.onBondClicked = function(bond) { ... };
editor.onAtomHighlight = function(atom) { ... };
editor.onBondHighlight = function(bond) { ... };
```

The full API is in [USAGE.md](USAGE.md). Every public function has
JSDoc comments in `editor/MolEditor.js`.

---

## 10. Capturing student work from a parent page

A common LMS / lecture-tool pattern: a student draws a molecule, the
parent page grabs the SMILES, scores it, and gives feedback.

```html
<!-- Parent page -->
<form id="quiz">
    <p>Q: Draw the structure of acetic acid.</p>
    <div id="quiz-editor" style="width:100%;height:400px"></div>
    <button type="button" id="check">Check answer</button>
    <p id="result"></p>
</form>

<script src="https://asad.github.io/bime-dist/dist/bime.min.js"></script>
<script>
    var ed = new MolEditor('quiz-editor', '100%', '400px');

    document.getElementById('check').onclick = function() {
        var studentSmiles = ed.getCanonicalSmiles();
        var expected = SmilesWriter.write(SmilesParser.parse('CC(=O)O'));
        if (studentSmiles === expected) {
            document.getElementById('result').textContent = '✓ Correct!';
        } else {
            document.getElementById('result').textContent = '✗ Try again.';
        }
    };
</script>
```

For partial credit / scaffold-aware grading, use SMARTS:

```js
var hasCarboxyl = SmartsMatch.hasMatch(
    ed.getMolecule(),
    SmartsParser.parse('[CX3](=O)O')
);
```

For LTI / xAPI integration with an LMS gradebook, capture the SMILES
on submit and POST it to your scoring endpoint. BIME itself is
client-only — it does not call your API; you do.

---

## 11. CSP and security considerations

BIME's default workbench ships with a strict Content Security Policy:

```
default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'
```

When embedding BIME in **your** page, **your** CSP applies. Make sure:

- `script-src` allows the origin you're loading BIME from (e.g.
  `'self'`, or `https://asad.github.io` if you cross-origin-load).
- `connect-src 'self'` is fine — BIME makes zero outbound network
  calls.
- If you use the bundle with SRI (§3), `script-src` must allow
  `crossorigin="anonymous"` and the integrity hash — which it does
  by default.

### XSS-safe DOM construction

BIME never uses `innerHTML` on user-controlled data. All atom labels
go through `textContent` or `createElementNS` for SVG. Halo highlights
use `createElementNS` SVG circles — strict-CSP-safe.

If you extend BIME with custom UI, keep the same discipline: never
splice user-supplied SMILES into `innerHTML`.

### `frame-ancestors` and iframe embedding

If you host BIME yourself and want to allow embedding via iframe (§2),
either:

1. Edit the `frame-ancestors` directive in the HTML `<meta>` to allow
   your LMS / blog domain, **or**
2. Set the `Content-Security-Policy` HTTP header at the server level
   to override the in-page meta.

The default `frame-ancestors 'none'` blocks all framing for security.
The official live demo at `asad.github.io/bime-dist/` uses this default
— so iframe embeds from third-party domains may need to use a
self-hosted copy with relaxed `frame-ancestors`.

---

## More resources

- [USAGE.md](USAGE.md) — the full BIME programmatic API.
- [HOSTING.md](HOSTING.md) — deploying BIME on your own infrastructure.
- [EDUCATORS.md](EDUCATORS.md) / [STUDENTS.md](STUDENTS.md) /
  [RESEARCHERS.md](RESEARCHERS.md) — audience-specific guides.
- [GitHub Issues](https://github.com/asad/bime-dist/issues) — questions,
  bug reports, feature requests.
