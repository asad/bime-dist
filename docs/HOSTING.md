# Self-hosting BIME

BIME is a **pure-JavaScript, zero-dependency** static-file application.
It runs entirely in the browser. There is no backend, no database, and no
server-side computation.

This means you can host BIME on **anything** that serves static files —
from a USB stick to a global CDN.

This guide covers the most common deployment options.

---

## Table of contents

1. [Just run it locally (no server)](#1-just-run-it-locally-no-server)
2. [Local web server (one command)](#2-local-web-server-one-command)
3. [GitHub Pages](#3-github-pages)
4. [Netlify](#4-netlify)
5. [Cloudflare Pages](#5-cloudflare-pages)
6. [Apache (university web space, shared hosting)](#6-apache-university-web-space-shared-hosting)
7. [Nginx](#7-nginx)
8. [Docker](#8-docker)
9. [USB stick / shared drive / intranet](#9-usb-stick--shared-drive--intranet)
10. [IPFS](#10-ipfs)
11. [Recommended HTTP headers (CSP, integrity, caching)](#11-recommended-http-headers-csp-integrity-caching)
12. [Customising what gets deployed](#12-customising-what-gets-deployed)
13. [Updating to a new BIME version](#13-updating-to-a-new-bime-version)

---

## 1. Just run it locally (no server)

The simplest option. BIME works directly from `file://` URIs.

```bash
git clone https://github.com/asad/bime-dist.git
cd bime-dist

# macOS
open workbench.html

# Linux
xdg-open workbench.html

# Windows
start workbench.html
```

**Caveats** (they almost never matter for normal use, but they exist):

- A few browsers (Safari in particular) restrict `fetch()` from
  `file://` URIs. BIME's default workbench does not use `fetch` at all,
  so this is fine. If you customise BIME and add `fetch`, switch to a
  local web server (next section).
- Service workers do not run from `file://`. BIME does not use service
  workers. Again, only relevant if you customise.

---

## 2. Local web server (one command)

Useful if you want to preview the site exactly as it would behave in
production, or if you customise BIME and want to use `fetch` / service
workers.

### Python (built into macOS, Linux, and most Windows installs)

```bash
cd bime-dist
python3 -m http.server 8080
# open http://localhost:8080/workbench.html
```

### Node (no install needed if Node.js is present)

```bash
npx serve -p 8080
# or
npx http-server -p 8080
```

### Ruby

```bash
ruby -run -e httpd . -p 8080
```

### PHP

```bash
php -S localhost:8080
```

Any of these will serve BIME with the default MIME types and let you
test as if it were on a real server.

---

## 3. GitHub Pages

This is how the official live demo at https://asad.github.io/bime-dist/
is hosted. It is free for public repositories and supports custom
domains.

### Steps

1. **Fork** [`asad/bime-dist`](https://github.com/asad/bime-dist) (or
   create your own repo with the BIME files).
2. Go to **Settings → Pages**.
3. Under **Source**, choose **Deploy from a branch**.
4. Pick **`main`** branch, folder **`/ (root)`**, and click **Save**.
5. Wait ~30 seconds; your fork is live at
   `https://<your-username>.github.io/<repo-name>/workbench.html`.

### Custom domain

1. Add a `CNAME` file at the repo root with one line: your domain
   (e.g. `chem.example.edu`).
2. In your DNS provider, add a `CNAME` record pointing your subdomain
   to `<your-username>.github.io`.
3. In **Settings → Pages**, enable **Enforce HTTPS**.

GitHub provides a free Let's Encrypt certificate. SSL takes 24-48 hours
the first time but is automatic thereafter.

### Important: keep `.nojekyll`

BIME repositories include a top-level `.nojekyll` file. **Do not delete
it** — it tells GitHub Pages to skip Jekyll processing, which would
otherwise mangle filenames starting with `_`.

---

## 4. Netlify

Drag-and-drop simplicity, free tier with HTTPS and a `*.netlify.app`
domain.

### One-click deploy

1. Go to https://app.netlify.com/drop.
2. Drag the `bime-dist` folder onto the page.
3. Done. Your site is live.

### Git-based deploy

1. Push the BIME repo to GitHub / GitLab / Bitbucket.
2. **New site → Import from Git → pick the repo.**
3. **Build command:** *(leave empty — no build needed)*.
4. **Publish directory:** `.` (the repo root).
5. Click **Deploy**.

### Custom domain

Add your domain under **Domain settings** and follow the DNS
instructions. Netlify provisions a Let's Encrypt cert automatically.

---

## 5. Cloudflare Pages

Free, fast, with a generous bandwidth allowance.

1. Push BIME to a Git provider.
2. **Cloudflare Pages → Create a project → Connect to Git.**
3. **Build command:** *(empty)*.
4. **Build output directory:** `/`.
5. Deploy. You get a `*.pages.dev` URL immediately.

For a custom domain, attach it under **Custom domains**.

---

## 6. Apache (university web space, shared hosting)

Many universities give staff and students static-file web space served
by Apache. BIME drops into any such space without configuration.

### Steps

1. Copy the entire BIME directory tree to your web space (typically
   `~/public_html/` or `~/www/`).
2. Browse to `https://<your-username>.<your-uni>.edu/~<your-username>/bime/workbench.html`.

### Optional: a virtual host

If you have full server access:

```apache
<VirtualHost *:443>
    ServerName chem.example.edu
    DocumentRoot /var/www/bime

    <Directory /var/www/bime>
        Options Indexes FollowSymLinks
        AllowOverride None
        Require all granted
    </Directory>

    # Recommended security headers (see §11 below)
    Header always set X-Content-Type-Options "nosniff"
    Header always set Referrer-Policy "strict-origin-when-cross-origin"
    Header always set X-Frame-Options "DENY"
    # CSP is already set by BIME's HTML <meta> tag; no need to override.

    # SSL configured via your existing Let's Encrypt / institutional cert
    SSLEngine on
    SSLCertificateFile /etc/letsencrypt/live/chem.example.edu/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/chem.example.edu/privkey.pem
</VirtualHost>
```

---

## 7. Nginx

```nginx
server {
    listen 443 ssl http2;
    server_name chem.example.edu;

    root /var/www/bime;
    index index.html;

    # Security headers
    add_header X-Content-Type-Options nosniff always;
    add_header Referrer-Policy "strict-origin-when-cross-origin" always;
    add_header X-Frame-Options DENY always;

    # Brotli / gzip for smaller transfers
    gzip on;
    gzip_types text/css application/javascript text/html application/json image/svg+xml;
    gzip_min_length 1024;

    # Long-cache for hashed bundles, short-cache for HTML
    location ~* \.(js|css|svg|png|jpg|woff2)$ {
        expires 30d;
        add_header Cache-Control "public, max-age=2592000, immutable";
    }
    location ~* \.html$ {
        expires 1h;
        add_header Cache-Control "public, max-age=3600";
    }

    location / {
        try_files $uri $uri/ =404;
    }

    ssl_certificate     /etc/letsencrypt/live/chem.example.edu/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/chem.example.edu/privkey.pem;
}

# HTTP -> HTTPS redirect
server {
    listen 80;
    server_name chem.example.edu;
    return 301 https://$server_name$request_uri;
}
```

---

## 8. Docker

A minimal Dockerfile using Nginx:

```dockerfile
FROM nginx:1.27-alpine
COPY . /usr/share/nginx/html
EXPOSE 80
```

Build and run:

```bash
docker build -t bime .
docker run -p 8080:80 bime
# open http://localhost:8080/workbench.html
```

For Docker Compose with HTTPS via a reverse proxy (Caddy or Traefik),
see your reverse-proxy's documentation — BIME has no special
requirements.

---

## 9. USB stick / shared drive / intranet

BIME's "drop anywhere" property makes it perfect for environments
without internet:

- **Computer labs** with no outbound internet — copy the `bime-dist`
  folder to a shared drive; students open `workbench.html` directly.
- **Field work** — copy to a USB stick; works from any laptop.
- **Air-gapped environments** — the same applies.

No external network calls are made by default. You can verify this in
the browser DevTools **Network** tab — the only requests are for the
local files.

---

## 10. IPFS

For long-term, decentralised hosting:

```bash
ipfs add -r bime-dist
# IPFS returns a CID like QmXyz...

# Then access via any IPFS gateway:
https://ipfs.io/ipfs/QmXyz...
https://cloudflare-ipfs.com/ipfs/QmXyz...
```

You can also pin the CID with services like Pinata or Filebase to
ensure availability.

---

## 11. Recommended HTTP headers (CSP, integrity, caching)

BIME ships with a strict
[Content Security Policy](https://developer.mozilla.org/en-US/docs/Web/HTTP/CSP)
in the HTML `<meta>` tag of each page. You generally do **not** need to
override it. If you do set CSP at the HTTP-header level, ensure it
matches or is stricter than the in-page meta:

```
Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'
```

Other recommended headers (covered in the Apache and Nginx examples
above):

- `X-Content-Type-Options: nosniff` — prevents MIME-type confusion.
- `Referrer-Policy: strict-origin-when-cross-origin` — privacy.
- `X-Frame-Options: DENY` — clickjacking protection (CSP
  `frame-ancestors 'none'` covers this for modern browsers).
- `Strict-Transport-Security: max-age=31536000; includeSubDomains` —
  only set this if you have HTTPS configured properly and you are sure
  you'll keep it.

### Subresource Integrity (SRI)

If you serve BIME from one domain but want to load the bundle from
another (e.g. your own CDN), use the SRI hashes in `dist/SRI.txt`:

```html
<script src="https://your-cdn/bime.min.js"
        integrity="sha384-<value-from-dist/SRI.txt>"
        crossorigin="anonymous"></script>
```

The browser will refuse to execute the bundle if the hash does not
match — protection against CDN compromise or transport tampering.

---

## 12. Customising what gets deployed

You don't need every file. The minimum to ship the workbench is:

```
workbench.html        # The main app
editor/*.js           # Source modules (or use dist/bime.min.js instead)
common-molecules.js   # Built-in molecule database
css/style.css         # Stylesheet
js/nav.js             # Top-level page script
images/               # SVG icons + logo
```

You can omit:

- `tests/` and `test*.js` (only needed for development).
- `tools/` (only needed for building / auditing).
- `paper/` (research drafts; intentionally excluded from the public dist).
- `examples.html`, `screenshots.html`, `docs.html`, `index.html` — only
  needed if you want a multi-page site rather than just the workbench.

For a minimal embedded deployment, see [docs/EMBED.md](EMBED.md).

---

## 13. Updating to a new BIME version

When a new BIME release is available:

```bash
cd bime-dist
git fetch origin
git pull origin main
# or download the new release tarball from
# https://github.com/asad/bime-dist/releases
```

Then:

- **GitHub Pages / Netlify / Cloudflare** — auto-deploy on push.
- **Apache / Nginx / Docker** — `rsync` or `git pull` to your server,
  then bump the cache headers if you cache aggressively.

To pin to a specific version, check out the tag:

```bash
git checkout v1.4.4
```

This guarantees byte-for-byte reproducibility and is recommended for
production deployments.

### Verifying the bundle integrity

After upgrading:

```bash
cd dist
shasum -a 256 -c MANIFEST.sha256
```

This checks every file against the manifest and reports tampering.

---

## Need more?

- For technical embedding (LMS, blog, Jupyter): see [EMBED.md](EMBED.md).
- For programmatic control of BIME from your own JavaScript: see
  [USAGE.md](USAGE.md).
- Open an issue if your hosting environment isn't covered:
  https://github.com/asad/bime-dist/issues
