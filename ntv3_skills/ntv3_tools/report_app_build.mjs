// Build the NTv3 report app (the Vite UI with GenomeBrowser) into one HTML
// shell with an empty report slot. `ntv3 report` fills the slot per run, so a
// run page needs no dev server and no Node at report time.
//
// Usage: node report_app_build.mjs <app-dir> <out-shell.html>
// <app-dir> is the report app (index.html, src/, node_modules/ with vite).
import { mkdtempSync, readFileSync, readdirSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { pathToFileURL } from "node:url";

const [appDir, outFile] = process.argv.slice(2);
if (!appDir || !outFile) {
  console.error("usage: node report_app_build.mjs <app-dir> <out-shell.html>");
  process.exit(2);
}
const DATA_TAG = '<script type="application/json" id="ntv3-report">null</script>';

const { build } = await import(
  pathToFileURL(join(appDir, "node_modules", "vite", "dist", "node", "index.js")).href
);
const outDir = mkdtempSync(join(tmpdir(), "ntv3-report-app-"));
try {
  await build({
    root: appDir,
    configFile: join(appDir, "vite.config.ts"),
    base: "./",
    // public/ holds a demo report.json; the shell carries no run of its own.
    publicDir: false,
    logLevel: "warn",
    build: {
      outDir,
      emptyOutDir: true,
      cssCodeSplit: false,
      modulePreload: false,
      assetsInlineLimit: Number.MAX_SAFE_INTEGER,
      // One page, one script: the lazy GenomeBrowser chunk is folded in.
      chunkSizeWarningLimit: 4096,
      rolldownOptions: { output: { codeSplitting: false } },
    },
  });
  const assets = readdirSync(join(outDir, "assets"));
  const scripts = assets.filter((name) => name.endsWith(".js"));
  const styles = assets.filter((name) => name.endsWith(".css"));
  if (scripts.length !== 1) {
    throw new Error(`Expected one JS bundle, got ${scripts.length}: ${scripts.join(", ")}`);
  }
  const js = readFileSync(join(outDir, "assets", scripts[0]), "utf8").replace(
    /<\/script/gi,
    "<\\/script",
  );
  const css = styles
    .map((name) => readFileSync(join(outDir, "assets", name), "utf8"))
    .join("\n");
  let html = readFileSync(join(outDir, "index.html"), "utf8")
    // External fonts would need the network; the page falls back to system fonts.
    .replace(/\s*<link[^>]*fonts\.(googleapis|gstatic)\.com[^>]*>/g, "")
    .replace(/\s*<link rel="stylesheet"[^>]*href="\.\/assets\/[^"]+\.css"[^>]*>/g, "")
    .replace(/\s*<script type="module"[^>]*src="\.\/assets\/[^"]+\.js"[^>]*><\/script>/g, "");
  html = html
    .replace("</head>", () => `<style>\n${css}\n</style>\n</head>`)
    .replace("</body>", () => `${DATA_TAG}\n<script type="module">\n${js}\n</script>\n</body>`);
  if (html.includes('src="./assets/') || html.includes('href="./assets/')) {
    throw new Error("A built asset is still linked, not inlined.");
  }
  writeFileSync(outFile, html);
  console.log(`Wrote ${outFile} (${Buffer.byteLength(html)} bytes)`);
} finally {
  rmSync(outDir, { recursive: true, force: true });
}
