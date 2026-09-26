/** The editor's entry point: mount the app, and let it load the module itself. */

import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { App } from "./App";
import "./styles/index.css";

const root = document.getElementById("root");
if (root === null) {
  throw new Error("index.html has no #root to mount into");
}

createRoot(root).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
