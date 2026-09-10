import { Tabs } from 'expo-router';
import { Platform } from 'react-native';
import { useSafeAreaInsets } from 'react-native-safe-area-context';

import { HapticTab } from '@/components/haptic-tab';
import { IconSymbol } from '@/components/ui/icon-symbol';
import { useColors } from '@/hooks/use-colors';

export default function TabLayout() {
  const colors = useColors('dark');
  const insets = useSafeAreaInsets();
  const bottomPadding = Platform.OS === 'web' ? 10 : Math.max(insets.bottom, 8);

  return (
    <Tabs screenOptions={{
      headerShown: false,
      tabBarActiveTintColor: colors.primary,
      tabBarInactiveTintColor: '#647891',
      tabBarButton: HapticTab,
      tabBarStyle: { paddingTop: 7, paddingBottom: bottomPadding, height: 58 + bottomPadding, backgroundColor: '#0B1728', borderTopColor: '#22334A', borderTopWidth: 1 },
      tabBarLabelStyle: { fontSize: 9, fontWeight: '700' },
    }}>
      <Tabs.Screen name="index" options={{ title: 'Overview', tabBarIcon: ({ color }) => <IconSymbol name="square.grid.2x2.fill" size={20} color={color} /> }} />
      <Tabs.Screen name="gaps" options={{ title: 'Gaps', tabBarIcon: ({ color }) => <IconSymbol name="exclamationmark.triangle.fill" size={20} color={color} /> }} />
      <Tabs.Screen name="engines" options={{ title: 'Engines', tabBarIcon: ({ color }) => <IconSymbol name="cpu.fill" size={20} color={color} /> }} />
      <Tabs.Screen name="audit" options={{ title: 'Audit', tabBarIcon: ({ color }) => <IconSymbol name="list.bullet.rectangle.portrait.fill" size={20} color={color} /> }} />
    </Tabs>
  );
}
