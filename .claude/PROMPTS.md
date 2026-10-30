# Claude Code Prompts

## Start Session
```
Read .claude/CLAUDE.md and .claude/context/[BRANCH].md
Summarize: feature purpose, current progress, next steps.
```

## End Session  
```
Update .claude/context/[BRANCH].md:
- What we accomplished
- Updated status
- Next steps
- Important context for next session
```

## New Branch
```
Create .claude/context/[NEW-BRANCH].md from the template.
Analyze git diff to populate key files section.
```
