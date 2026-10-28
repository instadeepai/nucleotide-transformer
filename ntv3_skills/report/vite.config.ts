import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      buffer: "buffer/",
    },
  },
  server: {
    port: 5174,
    strictPort: true,
    fs: {
      // The repo root, so the dev server can read results/<run>/report.json.
      allow: ["..", "../.."],
    },
  },
});
