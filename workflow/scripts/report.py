"""Render a self-contained, escaped report with actual outputs and provenance."""
import argparse
import html
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    for name in ("title", "summary", "provenance", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    summary = json.loads(Path(args.summary).read_text())
    provenance = json.loads(Path(args.provenance).read_text())
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    title = html.escape(args.title)
    rows = "".join(f"<tr><th>{html.escape(k)}</th><td>{html.escape(str(v))}</td></tr>"
                   for k, v in summary.items() if k != "artifacts")
    links = []
    for name in summary.get("artifacts", []):
        if Path(name).name != name:
            raise ValueError("Report artifact links must be local basenames")
        links.append(f'<li><a href="{html.escape(name, quote=True)}">{html.escape(name)}</a></li>')
    chart = ""
    figure = Path(args.summary).parent / "panel_counts.svg"
    if summary.get("synthetic") and figure.exists():
        # This SVG is produced by the local, reviewed demo script, not user HTML.
        chart = figure.read_text()
    doc = f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title><style>
body{{font:16px/1.6 system-ui,sans-serif;color:#18343d;background:#f5f7f6;margin:0}}
main{{max-width:1000px;margin:40px auto;padding:36px;background:white;border-radius:12px}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%}}
th,td{{padding:9px;border-bottom:1px solid #dce5e4;text-align:left;vertical-align:top}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f4f4;padding:18px}}
svg{{max-width:100%;height:auto}}a{{color:#006b83}}
</style><main><h1>{title}</h1><p>Generated from declared inputs by Snakemake.
Read the limitations and provenance before interpreting these outputs.</p>
{chart}<h2>Results</h2><table>{rows}</table><h2>Output files</h2><ul>{''.join(links)}</ul>
<h2>Parameters and provenance</h2><pre>{html.escape(json.dumps(provenance, indent=2, sort_keys=True))}</pre>
</main></html>'''
    out.write_text(doc)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
