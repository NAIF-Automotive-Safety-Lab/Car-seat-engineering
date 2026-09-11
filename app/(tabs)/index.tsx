import { useRouter } from 'expo-router';
import { Platform, Pressable, ScrollView, StyleSheet, Text, View } from 'react-native';
import * as Haptics from 'expo-haptics';

import { ScreenContainer } from '@/components/screen-container';
import { IconSymbol } from '@/components/ui/icon-symbol';
import { useColors } from '@/hooks/use-colors';

type Tone = 'mint' | 'amber' | 'coral' | 'slate';

const tones: Record<Tone, string> = {
  mint: '#2DE0B2',
  amber: '#FFBF69',
  coral: '#FF6B76',
  slate: '#8DA0B8',
};

type ControlDomain = { label: string; state: string; detail: string; tone: Tone; icon: 'shield.lefthalf.filled' | 'account-tree.fill' | 'cpu.fill' | 'doc.text.magnifyingglass' | 'list.bullet.rectangle.portrait.fill' | 'waveform.path.ecg' };
const controlDomains: ControlDomain[] = [
  { label: 'Trust kernel', state: 'GATED', detail: 'No execution without exact inputs', tone: 'mint', icon: 'shield.lefthalf.filled' },
  { label: 'Repository', state: 'BLOCKED', detail: 'Canonical identity required', tone: 'coral', icon: 'account-tree.fill' },
  { label: 'Engine adapters', state: 'UNKNOWN', detail: 'Adapter presence ≠ solver success', tone: 'amber', icon: 'cpu.fill' },
  { label: 'Evidence', state: 'GATED', detail: 'Provenance-bound outputs only', tone: 'amber', icon: 'doc.text.magnifyingglass' },
  { label: 'Audit manager', state: 'READY', detail: 'Every gate decision recorded', tone: 'mint', icon: 'list.bullet.rectangle.portrait.fill' },
  { label: 'Drift manager', state: 'NOT_PROVEN', detail: 'Baseline hash is unavailable', tone: 'slate', icon: 'waveform.path.ecg' },
];

const kernelModules = [
  ['Artifact manager', 'SHA and source ownership'],
  ['Test registry', 'Versioned execution contracts'],
  ['Execution bus', 'Queue only after trust gate'],
  ['Reproducibility', 'Compare inputs and outputs'],
] as const;

export default function HomeScreen() {
  const router = useRouter();
  const colors = useColors('dark');

  const pulse = async () => {
    if (Platform.OS !== 'web') await Haptics.impactAsync(Haptics.ImpactFeedbackStyle.Light);
  };

  return (
    <ScreenContainer containerClassName="bg-[#08111F]" className="px-5" edges={['top', 'left', 'right']}>
      <ScrollView contentContainerStyle={styles.scrollContent} showsVerticalScrollIndicator={false}>
        <View style={styles.headerRow}>
          <View style={styles.brandRow}>
            <View style={styles.brandMark}><IconSymbol name="shield.lefthalf.filled" size={21} color="#08111F" /></View>
            <View><Text style={styles.eyebrow}>V7-AEGIS / EXECUTION KERNEL</Text><Text style={styles.brandName}>AEGIS-X COMMAND CENTER</Text></View>
          </View>
          <Pressable onPress={pulse} style={({ pressed }) => [styles.iconButton, pressed && styles.pressed]}><IconSymbol name="arrow.clockwise" size={18} color={colors.foreground} /></Pressable>
        </View>

        <View style={styles.kernelBanner}>
          <View style={styles.bannerTop}><View style={styles.liveDot} /><Text style={styles.bannerLabel}>KERNEL ONLINE · TRUST GATED</Text><Text style={styles.bannerState}>NO FAKE PASS</Text></View>
          <Text style={styles.bannerTitle}>Single execution and qualification heart</Text>
          <Text style={styles.bannerCopy}>AEGIS-X plans, blocks, records, and audits every future CAD / CAE test. A registered adapter is not evidence that a solver succeeded.</Text>
          <View style={styles.bannerRule}><View style={styles.ruleSegment} /><Text style={styles.ruleText}>BLOCK BEFORE EXECUTION</Text></View>
        </View>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>Kernel health</Text><Text style={styles.sectionSubtitle}>Truth-bound operational domains</Text></View><Text style={styles.sectionStatus}>LIVE POLICY</Text></View>
        <View style={styles.healthGrid}>
          {controlDomains.map((item) => (
            <View key={item.label} style={styles.healthCard}>
              <View style={styles.healthTop}><IconSymbol name={item.icon} size={17} color={tones[item.tone]} /><View style={[styles.healthIndicator, { backgroundColor: tones[item.tone] }]} /></View>
              <Text style={styles.healthLabel}>{item.label}</Text>
              <Text style={[styles.healthState, { color: tones[item.tone] }]}>{item.state}</Text>
              <Text style={styles.healthDetail}>{item.detail}</Text>
            </View>
          ))}
        </View>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>R4.1 qualification</Text><Text style={styles.sectionSubtitle}>First-class workflow, not a viewer shortcut</Text></View><Text style={styles.blockedText}>BLOCKED</Text></View>
        <Pressable onPress={() => { void pulse(); router.push('/gaps' as never); }} style={({ pressed }) => [styles.qualificationCard, pressed && styles.pressed]}>
          <View style={styles.qualificationIcon}><IconSymbol name="exclamationmark.triangle.fill" size={20} color="#FFBF69" /></View>
          <View style={styles.qualificationBody}><View style={styles.qualificationTop}><Text style={styles.qualificationEyebrow}>R4.1-QUALIFICATION · v1.0.0</Text><Text style={styles.blockedText}>GATE</Text></View><Text style={styles.qualificationTitle}>Identity → SHA → STEP → B-Rep → Evidence</Text><Text style={styles.qualificationCopy}>Missing canonical artifact proof keeps every downstream check blocked. No computed result is promoted to physical validation.</Text><Text style={styles.linkText}>OPEN GAPS  ›</Text></View>
        </Pressable>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>Kernel modules</Text><Text style={styles.sectionSubtitle}>Centralized capabilities, one execution path</Text></View></View>
        <View style={styles.moduleGrid}>
          {kernelModules.map(([label, detail]) => <View key={label} style={styles.moduleCard}><View style={styles.moduleIcon}><IconSymbol name="chevron.right" size={15} color="#2DE0B2" /></View><Text style={styles.moduleLabel}>{label}</Text><Text style={styles.moduleDetail}>{detail}</Text></View>)}
        </View>

        <View style={styles.actionRow}>
          <Pressable onPress={() => { void pulse(); router.push('/engines' as never); }} style={({ pressed }) => [styles.actionCard, pressed && styles.pressed]}><IconSymbol name="cpu.fill" size={20} color="#2DE0B2" /><Text style={styles.actionTitle}>Engine registry</Text><Text style={styles.actionCopy}>OpenCascade · OpenRadioss · CalculiX · MBD adapters</Text></Pressable>
          <Pressable onPress={() => { void pulse(); router.push('/audit' as never); }} style={({ pressed }) => [styles.actionCard, pressed && styles.pressed]}><IconSymbol name="list.bullet.rectangle.portrait.fill" size={20} color="#FFBF69" /><Text style={styles.actionTitle}>Audit + evidence</Text><Text style={styles.actionCopy}>Trace provenance, gates, drift, and reproducibility</Text></Pressable>
        </View>

        <Text style={styles.footer}>ASSUMPTION ≠ MEASUREMENT · MODEL RESULT ≠ PHYSICAL VALIDATION · BLOCKED ≠ PASS</Text>
      </ScrollView>
    </ScreenContainer>
  );
}

const styles = StyleSheet.create({
  scrollContent: { paddingTop: 14, paddingBottom: 42, gap: 18 },
  headerRow: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'center' },
  brandRow: { flexDirection: 'row', alignItems: 'center', gap: 10 },
  brandMark: { width: 39, height: 39, borderRadius: 13, backgroundColor: '#2DE0B2', alignItems: 'center', justifyContent: 'center' },
  eyebrow: { color: '#8192A9', fontSize: 8, fontWeight: '800', letterSpacing: 1.25 },
  brandName: { color: '#F5F8FC', fontSize: 16, fontWeight: '900', letterSpacing: 0.7, marginTop: 3 },
  iconButton: { width: 40, height: 40, borderRadius: 13, backgroundColor: '#101D31', borderWidth: 1, borderColor: '#22334A', alignItems: 'center', justifyContent: 'center' },
  pressed: { opacity: 0.72, transform: [{ scale: 0.98 }] },
  kernelBanner: { backgroundColor: '#12243A', borderRadius: 21, borderWidth: 1, borderColor: '#2B5160', padding: 18, overflow: 'hidden' },
  bannerTop: { flexDirection: 'row', alignItems: 'center', gap: 7 },
  liveDot: { width: 7, height: 7, borderRadius: 4, backgroundColor: '#2DE0B2' },
  bannerLabel: { color: '#2DE0B2', fontSize: 9, fontWeight: '900', letterSpacing: 1.05, flex: 1 },
  bannerState: { color: '#FFBF69', fontSize: 8, fontWeight: '900', letterSpacing: 0.8 },
  bannerTitle: { color: '#F5F8FC', fontSize: 22, fontWeight: '900', letterSpacing: -0.5, marginTop: 17 },
  bannerCopy: { color: '#B0C0D3', fontSize: 12, lineHeight: 18, marginTop: 8 },
  bannerRule: { flexDirection: 'row', alignItems: 'center', gap: 9, marginTop: 17 },
  ruleSegment: { width: 30, height: 2, backgroundColor: '#FF6B76' },
  ruleText: { color: '#FF6B76', fontSize: 9, fontWeight: '900', letterSpacing: 1 },
  sectionHeader: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'flex-end', marginTop: 2 },
  sectionTitle: { color: '#F5F8FC', fontSize: 17, fontWeight: '900' },
  sectionSubtitle: { color: '#8192A9', fontSize: 10, marginTop: 4 },
  sectionStatus: { color: '#2DE0B2', fontSize: 9, fontWeight: '900', letterSpacing: 1 },
  blockedText: { color: '#FF6B76', fontSize: 9, fontWeight: '900', letterSpacing: 0.8 },
  healthGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 9 },
  healthCard: { width: '31.8%', minWidth: 95, backgroundColor: '#101D31', borderRadius: 15, borderWidth: 1, borderColor: '#22334A', padding: 11 },
  healthTop: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'center' },
  healthIndicator: { width: 5, height: 5, borderRadius: 3 },
  healthLabel: { color: '#9CAFC1', fontSize: 9, fontWeight: '700', marginTop: 13 },
  healthState: { fontSize: 11, fontWeight: '900', marginTop: 5 },
  healthDetail: { color: '#70829A', fontSize: 8, lineHeight: 12, marginTop: 4 },
  qualificationCard: { flexDirection: 'row', gap: 12, backgroundColor: '#171D2B', borderRadius: 17, borderWidth: 1, borderColor: '#44394A', padding: 14 },
  qualificationIcon: { width: 38, height: 38, borderRadius: 12, backgroundColor: '#332C26', alignItems: 'center', justifyContent: 'center' },
  qualificationBody: { flex: 1 },
  qualificationTop: { flexDirection: 'row', justifyContent: 'space-between', gap: 8 },
  qualificationEyebrow: { color: '#8D9FB5', fontSize: 8, fontWeight: '800', letterSpacing: 0.8 },
  qualificationTitle: { color: '#F5F8FC', fontSize: 14, fontWeight: '900', marginTop: 8 },
  qualificationCopy: { color: '#A2B1C2', fontSize: 10, lineHeight: 16, marginTop: 6 },
  linkText: { color: '#FFBF69', fontSize: 9, fontWeight: '900', letterSpacing: 0.8, marginTop: 10 },
  moduleGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 9 },
  moduleCard: { width: '48.6%', backgroundColor: '#0E1A2B', borderRadius: 15, borderWidth: 1, borderColor: '#22334A', padding: 12 },
  moduleIcon: { width: 27, height: 27, borderRadius: 9, backgroundColor: '#173833', alignItems: 'center', justifyContent: 'center', marginBottom: 10 },
  moduleLabel: { color: '#E4ECF4', fontSize: 11, fontWeight: '900' },
  moduleDetail: { color: '#778AA2', fontSize: 9, lineHeight: 13, marginTop: 4 },
  actionRow: { flexDirection: 'row', gap: 9 },
  actionCard: { flex: 1, backgroundColor: '#101D31', borderRadius: 16, borderWidth: 1, borderColor: '#22334A', padding: 13 },
  actionTitle: { color: '#F5F8FC', fontSize: 12, fontWeight: '900', marginTop: 11 },
  actionCopy: { color: '#8192A9', fontSize: 9, lineHeight: 14, marginTop: 4 },
  footer: { color: '#5D6E84', fontSize: 8, fontWeight: '800', letterSpacing: 0.7, textAlign: 'center', lineHeight: 13, marginTop: 2 },
});
