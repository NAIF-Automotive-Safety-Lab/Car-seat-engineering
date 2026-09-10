import AsyncStorage from '@react-native-async-storage/async-storage';
import { useRouter } from 'expo-router';
import { useEffect, useMemo, useState } from 'react';
import { Platform, Pressable, ScrollView, StyleSheet, Text, View } from 'react-native';
import * as Haptics from 'expo-haptics';

import { ScreenContainer } from '@/components/screen-container';
import { IconSymbol } from '@/components/ui/icon-symbol';
import { useColors } from '@/hooks/use-colors';

const toneColors: Record<'success' | 'warning' | 'error', string> = { success: '#2DE0B2', warning: '#FFBF69', error: '#FF6B76' };
const healthItems: Array<{ label: string; value: string; detail: string; icon: 'account-tree.fill' | 'arrow.triangle.2.circlepath' | 'waveform.path.ecg' | 'cpu.fill' | 'checkmark.seal.fill' | 'doc.text.magnifyingglass'; tone: keyof typeof toneColors }> = [
  { label: 'Repository', value: 'BLOCKED', detail: 'Identity unverified', icon: 'account-tree.fill', tone: 'error' },
  { label: 'Sync gate', value: 'BLOCKED', detail: 'No remote proof', icon: 'arrow.triangle.2.circlepath', tone: 'error' },
  { label: 'Drift gate', value: 'UNKNOWN', detail: 'Awaiting baseline', icon: 'waveform.path.ecg', tone: 'warning' },
  { label: 'Engines', value: '20 / 20', detail: 'Registered', icon: 'cpu.fill', tone: 'success' },
  { label: 'Test health', value: 'BLOCKED', detail: 'No verified run', icon: 'checkmark.seal.fill', tone: 'error' },
  { label: 'Evidence', value: 'GATED', detail: 'Source required', icon: 'doc.text.magnifyingglass', tone: 'warning' },
];

export default function HomeScreen() {
  const router = useRouter();
  const colors = useColors('dark');
  const [lastScan, setLastScan] = useState('Not run in this session');
  const [isScanning, setIsScanning] = useState(false);

  useEffect(() => { AsyncStorage.getItem('aegis.lastScan').then((value) => { if (value) setLastScan(value); }); }, []);
  const blockers = useMemo(() => healthItems.filter((item) => item.tone === 'error').length, []);

  const runScan = async () => {
    if (isScanning) return;
    setIsScanning(true);
    if (Platform.OS !== 'web') await Haptics.impactAsync(Haptics.ImpactFeedbackStyle.Light);
    setTimeout(async () => {
      const stamp = new Date().toISOString().replace('T', ' ').slice(0, 19) + 'Z';
      await AsyncStorage.setItem('aegis.lastScan', stamp);
      setLastScan(stamp);
      setIsScanning(false);
      if (Platform.OS !== 'web') await Haptics.notificationAsync(Haptics.NotificationFeedbackType.Warning);
    }, 650);
  };

  return (
    <ScreenContainer containerClassName="bg-[#08111F]" className="px-5" edges={['top', 'left', 'right']}>
      <ScrollView contentContainerStyle={styles.scrollContent} showsVerticalScrollIndicator={false}>
        <View style={styles.headerRow}>
          <View style={styles.brandRow}><View style={styles.brandMark}><IconSymbol name="shield.lefthalf.filled" size={20} color="#08111F" /></View><View><Text style={styles.eyebrow}>ENGINEERING CONTROL</Text><Text style={styles.brandName}>V7-AEGIS</Text></View></View>
          <Pressable onPress={runScan} style={({ pressed }) => [styles.iconButton, pressed && styles.pressed]}><IconSymbol name="arrow.clockwise" size={19} color={colors.foreground} /></Pressable>
        </View>
        <View style={styles.statusRow}><View style={styles.statusDot} /><Text style={styles.statusText}>{isScanning ? 'SELF-DIAGNOSTIC RUNNING' : 'SYSTEM READY WITH BLOCKERS'}</Text><Text style={styles.statusTime}>{blockers} blockers</Text></View>

        <View style={styles.heroCard}><View style={styles.heroGlow} /><View style={styles.heroTopLine}><Text style={styles.sectionKicker}>QUALIFICATION GATE</Text><View style={styles.blockedPill}><Text style={styles.blockedPillText}>BLOCKED</Text></View></View><Text style={styles.heroTitle}>R4.1 artifact qualification</Text><Text style={styles.heroCopy}>Execution is held until the canonical repository identity and artifact evidence are verified.</Text><View style={styles.heroDivider} /><View style={styles.heroMetaRow}><View><Text style={styles.metaLabel}>BASELINE</Text><Text style={styles.metaValue}>R4.1 · IMMUTABLE</Text></View><View><Text style={styles.metaLabel}>EVIDENCE</Text><Text style={styles.metaValue}>UNVERIFIED</Text></View></View></View>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>System health</Text><Text style={styles.sectionSubtitle}>Gates are truth-bound, not optimistic UI.</Text></View><Text style={styles.scanTime}>{lastScan === 'Not run in this session' ? 'LIVE' : 'SCANNED'}</Text></View>
        <View style={styles.healthGrid}>{healthItems.map((item) => <View key={item.label} style={styles.healthCard}><View style={styles.healthTop}><IconSymbol name={item.icon} size={17} color={toneColors[item.tone]} /><View style={[styles.healthIndicator, { backgroundColor: toneColors[item.tone] }]} /></View><Text style={styles.healthLabel}>{item.label}</Text><Text style={[styles.healthValue, { color: toneColors[item.tone] }]}>{item.value}</Text><Text style={styles.healthDetail}>{item.detail}</Text></View>)}</View>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>Priority attention</Text><Text style={styles.sectionSubtitle}>Resolve evidence blockers before remediation.</Text></View></View>
        <Pressable onPress={() => router.push('/gaps' as never)} style={({ pressed }) => [styles.priorityCard, pressed && styles.pressed]}><View style={styles.priorityIcon}><IconSymbol name="exclamationmark.triangle.fill" size={19} color="#FFBF69" /></View><View style={styles.priorityBody}><View style={styles.priorityTitleRow}><Text style={styles.priorityTitle}>GAP-0001 · repository identity</Text><Text style={styles.prioritySeverity}>P0</Text></View><Text style={styles.priorityCopy}>Exact branch, commit, and manifest hash are required before a qualification run.</Text><Text style={styles.priorityLink}>VIEW GAP CENTER  ›</Text></View></Pressable>

        <View style={styles.sectionHeader}><View><Text style={styles.sectionTitle}>Operator actions</Text><Text style={styles.sectionSubtitle}>Safe entry points for the next verified step.</Text></View></View>
        <View style={styles.actionRow}><Pressable onPress={() => router.push('/engines' as never)} style={({ pressed }) => [styles.actionCard, pressed && styles.pressed]}><View style={styles.actionIcon}><IconSymbol name="cpu.fill" size={20} color="#2DE0B2" /></View><Text style={styles.actionTitle}>Engine center</Text><Text style={styles.actionCopy}>20 registered engines</Text></Pressable><Pressable onPress={() => router.push('/audit' as never)} style={({ pressed }) => [styles.actionCard, pressed && styles.pressed]}><View style={[styles.actionIcon, { backgroundColor: '#2B2331' }]}><IconSymbol name="list.bullet.rectangle.portrait.fill" size={20} color="#FFBF69" /></View><Text style={styles.actionTitle}>Audit trail</Text><Text style={styles.actionCopy}>Trace every gate decision</Text></Pressable></View>
        <Text style={styles.footer}>LAST DIAGNOSTIC · {lastScan}</Text>
      </ScrollView>
    </ScreenContainer>
  );
}

const styles = StyleSheet.create({
  scrollContent: { paddingTop: 14, paddingBottom: 40, gap: 18 }, headerRow: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'center' }, brandRow: { flexDirection: 'row', alignItems: 'center', gap: 10 }, brandMark: { width: 38, height: 38, borderRadius: 12, backgroundColor: '#2DE0B2', alignItems: 'center', justifyContent: 'center' }, eyebrow: { color: '#8192A9', fontSize: 9, fontWeight: '800', letterSpacing: 1.5 }, brandName: { color: '#F5F8FC', fontSize: 18, fontWeight: '800', letterSpacing: 1.1, marginTop: 2 }, iconButton: { width: 40, height: 40, borderRadius: 13, backgroundColor: '#101D31', borderWidth: 1, borderColor: '#22334A', alignItems: 'center', justifyContent: 'center' }, pressed: { opacity: 0.72, transform: [{ scale: 0.98 }] }, statusRow: { flexDirection: 'row', alignItems: 'center', gap: 8, marginTop: 2 }, statusDot: { width: 7, height: 7, borderRadius: 4, backgroundColor: '#FFBF69' }, statusText: { color: '#FFBF69', fontSize: 10, fontWeight: '800', letterSpacing: 0.8, flex: 1 }, statusTime: { color: '#8192A9', fontSize: 10, fontWeight: '700' }, heroCard: { overflow: 'hidden', backgroundColor: '#12243A', borderRadius: 22, padding: 20, borderWidth: 1, borderColor: '#26405B' }, heroGlow: { position: 'absolute', width: 180, height: 180, borderRadius: 90, backgroundColor: '#153A48', right: -60, top: -80, opacity: 0.7 }, heroTopLine: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'center' }, sectionKicker: { color: '#7F9AB5', fontSize: 10, fontWeight: '800', letterSpacing: 1.4 }, blockedPill: { backgroundColor: '#412733', borderRadius: 8, paddingHorizontal: 9, paddingVertical: 5 }, blockedPillText: { color: '#FF6B76', fontSize: 10, fontWeight: '900', letterSpacing: 0.8 }, heroTitle: { color: '#F5F8FC', fontSize: 23, fontWeight: '800', marginTop: 18, letterSpacing: -0.4 }, heroCopy: { color: '#AFC0D4', fontSize: 13, lineHeight: 20, marginTop: 8, maxWidth: 320 }, heroDivider: { height: 1, backgroundColor: '#29445E', marginTop: 19, marginBottom: 14 }, heroMetaRow: { flexDirection: 'row', gap: 28 }, metaLabel: { color: '#7890A9', fontSize: 9, fontWeight: '800', letterSpacing: 1.2 }, metaValue: { color: '#DCE7F2', fontSize: 11, fontWeight: '700', marginTop: 5 }, sectionHeader: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'flex-end', marginTop: 4 }, sectionTitle: { color: '#F5F8FC', fontSize: 17, fontWeight: '800' }, sectionSubtitle: { color: '#8192A9', fontSize: 11, marginTop: 4 }, scanTime: { color: '#2DE0B2', fontSize: 10, fontWeight: '800', letterSpacing: 1 }, healthGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 10 }, healthCard: { width: '31.8%', minWidth: 95, backgroundColor: '#101D31', borderRadius: 15, borderWidth: 1, borderColor: '#22334A', padding: 12 }, healthTop: { flexDirection: 'row', justifyContent: 'space-between', alignItems: 'center' }, healthIndicator: { width: 5, height: 5, borderRadius: 3 }, healthLabel: { color: '#9CAFC1', fontSize: 10, fontWeight: '700', marginTop: 14 }, healthValue: { fontSize: 12, fontWeight: '900', marginTop: 5, letterSpacing: 0.2 }, healthDetail: { color: '#70829A', fontSize: 9, marginTop: 4, lineHeight: 13 }, priorityCard: { flexDirection: 'row', gap: 13, backgroundColor: '#171D2B', borderRadius: 17, borderWidth: 1, borderColor: '#44394A', padding: 15 }, priorityIcon: { width: 38, height: 38, borderRadius: 12, backgroundColor: '#332C26', alignItems: 'center', justifyContent: 'center' }, priorityBody: { flex: 1 }, priorityTitleRow: { flexDirection: 'row', justifyContent: 'space-between', gap: 8 }, priorityTitle: { flex: 1, color: '#F5F8FC', fontSize: 13, fontWeight: '800' }, prioritySeverity: { color: '#FF6B76', fontSize: 11, fontWeight: '900' }, priorityCopy: { color: '#9DAFC2', fontSize: 11, lineHeight: 17, marginTop: 6 }, priorityLink: { color: '#FFBF69', fontSize: 10, fontWeight: '900', letterSpacing: 0.7, marginTop: 10 }, actionRow: { flexDirection: 'row', gap: 10 }, actionCard: { flex: 1, backgroundColor: '#101D31', borderRadius: 16, borderWidth: 1, borderColor: '#22334A', padding: 14 }, actionIcon: { width: 37, height: 37, borderRadius: 11, backgroundColor: '#173833', alignItems: 'center', justifyContent: 'center', marginBottom: 12 }, actionTitle: { color: '#F5F8FC', fontSize: 13, fontWeight: '800' }, actionCopy: { color: '#8192A9', fontSize: 10, lineHeight: 15, marginTop: 4 }, footer: { textAlign: 'center', color: '#5D6E84', fontSize: 9, fontWeight: '700', letterSpacing: 0.9, marginTop: 4 },
});
