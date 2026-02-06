# Quick Start: Figure Commentary

## Current Status ✓

Your commentary system is now **fault-tolerant** and will work reliably.

## Quick Setup (Choose One)

### Option 1: Use Deterministic Commentary (Works Now)

**No setup required!** Your current config already works:

```yaml
commentary:
  enabled: yes
  backend: none
```

Run your pipeline normally:
```r
targets::tar_make()
```

### Option 2: Enable AI Commentary with Claude

1. **Install Python package:**
   ```bash
   pip install anthropic
   ```

2. **Get API key** from https://console.anthropic.com/

3. **Set API key:**
   ```r
   # Option A: In R session
   Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-your-key-here")

   # Option B: In .Renviron file (persistent)
   usethis::edit_r_environ()
   # Add line: ANTHROPIC_API_KEY=sk-ant-your-key-here
   ```

4. **Update config.yml:**
   ```yaml
   commentary:
     enabled: yes
     backend: claude  # ← Change from "none" to "claude"
     claude_model: claude-sonnet-4.5-20250514
   ```

5. **Run pipeline:**
   ```r
   targets::tar_make()
   ```

### Option 3: Enable AI Commentary with OpenAI

1. **Install Python package:**
   ```bash
   pip install openai
   ```

2. **Get API key** from https://platform.openai.com/

3. **Set API key:**
   ```r
   Sys.setenv(OPENAI_API_KEY = "sk-proj-your-key-here")
   ```

4. **Update config.yml:**
   ```yaml
   commentary:
     enabled: yes
     backend: openai  # ← Change from "none" to "openai"
     openai_model: gpt-4o
   ```

5. **Run pipeline:**
   ```r
   targets::tar_make()
   ```

## What Changed?

The system now:
- ✓ Always falls back to deterministic commentary if AI fails
- ✓ Shows clear warnings about missing API keys or packages
- ✓ Never crashes your pipeline due to commentary issues
- ✓ Provides detailed error messages for debugging

## Troubleshooting

### "anthropic package not installed"
```bash
pip install anthropic
```

### "ANTHROPIC_API_KEY not set"
```r
Sys.setenv(ANTHROPIC_API_KEY = "your-key")
```

### "Commentary generation failed"
**Don't worry!** The system automatically fell back to deterministic commentary. Your pipeline continues successfully.

## Testing

```r
# Test just the commentary target
targets::tar_make(commentary_tbl)

# Check the output
targets::tar_read(commentary_tbl)

# View generated commentary files
list.files("outputs/commentary", pattern = "\\.json$")
```

## Need Help?

See detailed audit report: `COMMENTARY_AUDIT_FIXES.md`
