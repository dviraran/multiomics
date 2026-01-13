#!/usr/bin/env python3
"""
Figure Commentary Generator using OpenAI GPT-4 Vision API

This script analyzes multi-omics integration pipeline figures using OpenAI's vision
capabilities and generates structured commentary including observations, potential
issues, and recommended next steps.

Usage:
    python figure_commentary_openai.py \
        --image path/to/figure.png \
        --figure_id mofa_variance_heatmap \
        --context_json path/to/context.json \
        --out_json path/to/output.json \
        --model gpt-4o \
        --max_tokens 1500

Environment Variables:
    OPENAI_API_KEY: Your OpenAI API key (required)
"""

import argparse
import base64
import json
import os
import sys
import time
from pathlib import Path

try:
    from openai import OpenAI
except ImportError:
    print("Error: openai package not installed. Run: pip install openai", file=sys.stderr)
    sys.exit(1)


# System prompt for multi-omics figure commentary
SYSTEM_PROMPT = """You are an expert bioinformatics analyst specializing in multi-omics data integration.
Your task is to analyze visualization figures from a multi-omics integration pipeline and provide
structured, accurate commentary.

CRITICAL INSTRUCTIONS:
1. ONLY describe patterns, trends, and features that you can DIRECTLY SEE in the figure
2. DO NOT invent specific numeric values, gene names, or sample names unless they are clearly visible
3. If something is unclear, unreadable, or ambiguous, explicitly say so
4. Be factual and objective - avoid over-interpretation
5. Focus on patterns that are actionable for biological interpretation
6. Consider the multi-omics context: relationships between transcriptomics, proteomics, and metabolomics

Multi-omics specific considerations:
- MOFA2: Factors capturing variance from multiple views indicate shared biological signal
- DIABLO: Supervised method - good separation may indicate overfitting without CV validation
- SNF: Unsupervised clustering may reveal novel patient subgroups
- Concordance: RNA-protein correlation ~0.3-0.6 is typical; lower values suggest post-transcriptional regulation

Your analysis should help researchers:
- Understand what the figure shows
- Identify potential quality issues or concerns
- Make informed decisions about multi-omics integration results

Always output valid JSON matching the required schema."""


def build_user_prompt(figure_id: str, context: dict) -> str:
    """Build the user prompt with figure context."""

    plot_type = context.get("plot_type", "unknown")
    title = context.get("title", figure_id)
    description = context.get("description", "")
    x_axis = context.get("x_axis", "")
    y_axis = context.get("y_axis", "")
    section = context.get("section", "")

    # Experiment context
    organism = context.get("organism", "unknown")
    n_samples = context.get("n_samples", "unknown")
    omics_present = context.get("omics_present", "")
    condition_col = context.get("condition_col", "")

    # Method-specific context
    extra_context = []
    if context.get("mofa_total_variance"):
        extra_context.append(f"MOFA total variance explained: {context['mofa_total_variance']}%")
    if context.get("diablo_cv_error") is not None:
        extra_context.append(f"DIABLO CV error rate: {context['diablo_cv_error']*100:.1f}%")
    if context.get("snf_nmi") is not None:
        extra_context.append(f"SNF NMI with condition: {context['snf_nmi']:.3f}")
    if context.get("concordance_mean_cor") is not None:
        extra_context.append(f"Mean RNA-protein correlation: {context['concordance_mean_cor']:.3f}")

    prompt = f"""Analyze this figure from a multi-omics integration pipeline.

FIGURE CONTEXT:
- Figure ID: {figure_id}
- Plot Type: {plot_type}
- Section: {section}
- Title: {title}
- Description: {description}
- X-axis: {x_axis}
- Y-axis: {y_axis}

EXPERIMENT CONTEXT:
- Organism: {organism}
- Number of samples: {n_samples}
- Omics layers: {omics_present}
{f'- Condition variable: {condition_col}' if condition_col else ''}
{chr(10).join(f'- {ctx}' for ctx in extra_context) if extra_context else ''}

TASK:
Analyze the figure and provide commentary in the following JSON format:

{{
    "figure_id": "{figure_id}",
    "title": "Brief descriptive title",
    "what_is_this": "1-3 sentence description of what this figure shows and how to interpret it",
    "observations": [
        "Specific observation 1 based on visible patterns",
        "Specific observation 2",
        "... (3-8 observations)"
    ],
    "issues_checks": [
        "Potential issue or quality check 1",
        "... (0-5 items, empty list if none)"
    ],
    "next_steps": [
        "Recommended action 1 based on observations",
        "... (0-5 items)"
    ],
    "confidence": "high|medium|low",
    "limitations": "Any limitations in the analysis or things that were unclear"
}}

IMPORTANT REMINDERS:
- Only describe what you can actually see in the figure
- If labels are unreadable, say so rather than guessing
- Do not make up sample names, gene names, or numeric values
- Focus on patterns (clustering, separation, correlations, outliers)
- For integration plots: consider what the pattern means for multi-omics biology
- Be specific but avoid over-interpretation

Return ONLY the JSON object, no additional text."""

    return prompt


def encode_image(image_path: str) -> tuple[str, str]:
    """Encode image to base64 and determine media type."""
    path = Path(image_path)
    suffix = path.suffix.lower()

    media_types = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".gif": "image/gif",
        ".webp": "image/webp"
    }

    media_type = media_types.get(suffix, "image/png")

    with open(image_path, "rb") as f:
        image_data = base64.standard_b64encode(f.read()).decode("utf-8")

    return image_data, media_type


def call_openai_vision(
    image_path: str,
    figure_id: str,
    context: dict,
    model: str = "gpt-4o",
    max_tokens: int = 1500,
    max_retries: int = 3,
    retry_delay: float = 2.0
) -> dict:
    """Call OpenAI Vision API with retry logic."""

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("OPENAI_API_KEY environment variable not set")

    client = OpenAI(api_key=api_key)

    image_data, media_type = encode_image(image_path)
    user_prompt = build_user_prompt(figure_id, context)

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {
            "role": "user",
            "content": [
                {
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:{media_type};base64,{image_data}",
                        "detail": "high"
                    }
                },
                {
                    "type": "text",
                    "text": user_prompt
                }
            ]
        }
    ]

    last_error = None
    for attempt in range(max_retries):
        try:
            response = client.chat.completions.create(
                model=model,
                max_tokens=max_tokens,
                messages=messages
            )

            response_text = response.choices[0].message.content

            # Parse JSON
            if "```json" in response_text:
                json_start = response_text.find("```json") + 7
                json_end = response_text.find("```", json_start)
                response_text = response_text[json_start:json_end].strip()
            elif "```" in response_text:
                json_start = response_text.find("```") + 3
                json_end = response_text.find("```", json_start)
                response_text = response_text[json_start:json_end].strip()

            commentary = json.loads(response_text)

            # Validate required fields
            required_fields = ["figure_id", "what_is_this", "observations"]
            for field in required_fields:
                if field not in commentary:
                    commentary[field] = "" if field != "observations" else []

            for field in ["observations", "issues_checks", "next_steps"]:
                if field in commentary and not isinstance(commentary[field], list):
                    commentary[field] = [commentary[field]] if commentary[field] else []

            return commentary

        except Exception as e:
            last_error = e
            error_str = str(e).lower()

            if "rate" in error_str or "limit" in error_str or "429" in error_str:
                wait_time = retry_delay * (2 ** attempt)
                print(f"Rate limited. Waiting {wait_time}s before retry {attempt + 1}/{max_retries}", file=sys.stderr)
                time.sleep(wait_time)
            elif attempt < max_retries - 1:
                print(f"API error: {e}. Retrying {attempt + 1}/{max_retries}", file=sys.stderr)
                time.sleep(retry_delay)

    if isinstance(last_error, json.JSONDecodeError):
        return {
            "figure_id": figure_id,
            "title": "Commentary Parse Error",
            "what_is_this": "Unable to parse the AI-generated commentary.",
            "observations": ["Raw response could not be parsed as JSON"],
            "issues_checks": [],
            "next_steps": [],
            "confidence": "none",
            "limitations": f"JSON parse error: {str(last_error)}"
        }

    raise RuntimeError(f"Failed after {max_retries} attempts. Last error: {last_error}")


def main():
    parser = argparse.ArgumentParser(description="Generate figure commentary using OpenAI GPT-4 Vision API")
    parser.add_argument("--image", required=True, help="Path to the figure image file")
    parser.add_argument("--figure_id", required=True, help="Unique identifier for the figure")
    parser.add_argument("--context_json", required=True, help="Path to JSON file with figure context")
    parser.add_argument("--out_json", required=True, help="Path for output JSON file")
    parser.add_argument("--model", default="gpt-4o", help="OpenAI model to use")
    parser.add_argument("--max_tokens", type=int, default=1500, help="Maximum tokens in response")
    parser.add_argument("--max_retries", type=int, default=3, help="Maximum retry attempts")

    args = parser.parse_args()

    if not os.path.exists(args.image):
        print(f"Error: Image file not found: {args.image}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.context_json):
        print(f"Error: Context file not found: {args.context_json}", file=sys.stderr)
        sys.exit(1)

    with open(args.context_json, "r") as f:
        context = json.load(f)

    try:
        commentary = call_openai_vision(
            image_path=args.image,
            figure_id=args.figure_id,
            context=context,
            model=args.model,
            max_tokens=args.max_tokens,
            max_retries=args.max_retries
        )

        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(args.out_json, "w") as f:
            json.dump(commentary, f, indent=2)

        print(f"Commentary saved to: {args.out_json}")

    except Exception as e:
        print(f"Error generating commentary: {e}", file=sys.stderr)

        error_commentary = {
            "figure_id": args.figure_id,
            "title": "Commentary Generation Failed",
            "what_is_this": "Unable to generate automated commentary.",
            "observations": [],
            "issues_checks": [],
            "next_steps": [],
            "confidence": "none",
            "limitations": str(e)
        }

        with open(args.out_json, "w") as f:
            json.dump(error_commentary, f, indent=2)

        sys.exit(1)


if __name__ == "__main__":
    main()
