go_file = '/Users/yeojingi/Documents/Github/m2g/examples/mmseqs2go/database/uniprot_sprot_func_gog'

go_info = {}
with open(go_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split(',')
        if len(parts) < 3:
            continue
        go_id = parts[0]
        description = parts[1]
        parent_info = parts[2] if len(parts) > 2 else ""
        category = parts[3] if len(parts) > 3 else "unknown"
        
        parents = []
        if parent_info:
            parent_ids = [p.strip() for p in parent_info.split(';')]
            for p_id in parent_ids:
                if p_id.startswith('GO:'):
                    parents.append(p_id)
        
        go_info[go_id] = {
            'name': description,
            'parents': parents,
            'category': category
        }

# GO:0005615와 GO:0005576의 전체 경로 추적
for target in ["GO:0005615", "GO:0005576"]:
    print(f"\n📍 Full path for {target}:")
    current = target
    for step in range(15):
        if current not in go_info:
            print(f"  ❌ Not found: {current}")
            break
        
        go_data = go_info[current]
        print(f"  {step}: {current} -> {go_data['name'][:50]} (parent: {go_data['parents']})")
        
        parents = go_data.get('parents', [])
        if not parents:
            print(f"  ✓ REACHED ROOT (no parents)")
            break
        
        current = parents[0]
