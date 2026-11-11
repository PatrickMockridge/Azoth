# Bundled faces

The editor's three layouts are deliberate evocations, and typography is most of the evocation.
The faces below are vendored as WOFF2 so the look is identical on Linux, macOS and Windows and
does not depend on a network. Licences are in `licenses/`; the root `NOTICE` is the attribution
record.

| Family | Files | Weight | Source | Licence | Changed |
| --- | --- | --- | --- | --- | --- |
| Sysfont | `sysfont.woff2` | 400 | https://madebyarne.github.io/sysfont/ (v1.0, `sysfont.zip`) | SIL OFL 1.1 | No |
| Wine Tahoma | `wine-tahoma.ttf`, `wine-tahoma-bold.ttf` | 400, 700 | https://github.com/wine-mirror/wine/tree/master/fonts (`tahoma.ttf`, `tahomabd.ttf`; commit at fetch time on `master`) | LGPL-2.1-or-later (Copyright 2004 Larry Snyder, based on Bitstream Vera Sans, Copyright 2003 Bitstream Inc.) | No |
| Rajdhani | `rajdhani-500.woff2`, `rajdhani-600.woff2` | 500, 600 | Google Fonts, `Rajdhani:v17` Latin subset (https://fonts.google.com/specimen/Rajdhani) | SIL OFL 1.1 | No |
| Share Tech Mono | `share-tech-mono.woff2` | 400 | Google Fonts, `Share+Tech+Mono:v16` Latin subset (https://fonts.google.com/specimen/Share+Tech+Mono) | SIL OFL 1.1 | No |
| IBM Plex Mono | `ibm-plex-mono-400.woff2` | 400 | Google Fonts, `IBM+Plex+Mono:v20` Latin subset (https://fonts.google.com/specimen/IBM+Plex+Mono) | SIL OFL 1.1 | No |

## Why Wine Tahoma stays TTF

Wine ships TTF, and the browser serves it directly. Converting it to WOFF2 would only change the
container, but the LGPL notes that a changed copy must still be sourceable; keeping the upstream
files byte-for-byte makes the licence question disappear. The face is only fetched when a
Windows XP layout is active.

## Why Wine Tahoma

Wine's Tahoma is the open metric-compatible replacement for the XP UI face, and it is the only
faithful option with a licence that permits redistribution. It is LGPL-2.1-or-later; the font is
unmodified apart from the container conversion, and the licence text and source URL travel with
it. If a future policy refuses LGPL assets, `Liberation Sans` (OFL) is the fallback, at the cost
of the XP face's narrower metrics.
