import { COOKIE_NAME } from "../shared/const.js";
import { getSessionCookieOptions } from "./_core/cookies";
import { systemRouter } from "./_core/systemRouter";
import { protectedProcedure, publicProcedure, router } from "./_core/trpc";
import { completeQualificationRun, createQualificationRun, getQualificationRun, listAudits } from "./qualification";
import { completeKernelRun, getKernelRun, listKernelStatus, planKernelRun, recordKernelDrift, recordKernelEvidence, registerKernelEngine, registerKernelTest } from "./kernel";
import { z } from "zod";

export const appRouter = router({
  // if you need to use socket.io, read and register route in server/_core/index.ts, all api should start with '/api/' so that the gateway can route correctly
  system: systemRouter,
  auth: router({
    me: publicProcedure.query((opts) => opts.ctx.user),
    logout: publicProcedure.mutation(({ ctx }) => {
      const cookieOptions = getSessionCookieOptions(ctx.req);
      ctx.res.clearCookie(COOKIE_NAME, { ...cookieOptions, maxAge: -1 });
      return {
        success: true,
      } as const;
    }),
  }),

  qualification: router({
    create: protectedProcedure
      .input(z.object({
        executor: z.string().min(1),
        testId: z.string().min(1),
        testVersion: z.string().min(1),
        repository: z.string().url(),
        branch: z.string().min(1),
        commit: z.string().min(7),
        modelId: z.string().optional(),
        modelVersion: z.string().optional(),
        modelSha256: z.string().optional(),
        inputs: z.unknown(),
        engine: z.string().min(1),
        engineVersion: z.string().min(1),
        command: z.string().min(1),
      }))
      .mutation(({ ctx, input }) => createQualificationRun({ userId: ctx.user.id, ...input })),
    complete: protectedProcedure
      .input(z.object({
        runId: z.string().min(1),
        result: z.unknown(),
        output: z.unknown(),
        outputHash: z.string().min(1),
        resultStatus: z.enum(["PASS", "FAIL", "BLOCKED", "NOT_PROVEN"]),
        evidenceIds: z.array(z.string()),
      }))
      .mutation(({ ctx, input }) => completeQualificationRun({ userId: ctx.user.id, ...input })),
    get: protectedProcedure
      .input(z.object({ runId: z.string().min(1) }))
      .query(({ ctx, input }) => getQualificationRun(ctx.user.id, input.runId)),
    audits: protectedProcedure
      .input(z.object({ runId: z.string().min(1).optional() }).optional())
      .query(({ ctx, input }) => listAudits(ctx.user.id, input?.runId)),
  }),

  kernel: router({
    status: protectedProcedure.query(() => listKernelStatus()),
    plan: protectedProcedure
      .input(z.object({
        testId: z.string().min(1),
        testVersion: z.string().min(1),
        repository: z.string().optional(),
        commit: z.string().optional(),
        artifactSha256: z.string().optional(),
        model: z.string().optional(),
        inputHashes: z.record(z.string(), z.string()).optional(),
        engineId: z.string().min(1),
        engineVersion: z.string().min(1),
      }))
      .mutation(({ ctx, input }) => planKernelRun({ userId: ctx.user.id, ...input })),
    getRun: protectedProcedure
      .input(z.object({ runId: z.string().min(1) }))
      .query(({ ctx, input }) => getKernelRun(ctx.user.id, input.runId)),
    complete: protectedProcedure
      .input(z.object({ runId: z.string().min(1), outputHash: z.string().min(1), result: z.unknown(), resultStatus: z.enum(["PASS", "FAIL", "BLOCKED", "NOT_PROVEN"]), evidenceIds: z.array(z.string()) }))
      .mutation(({ ctx, input }) => completeKernelRun({ userId: ctx.user.id, ...input })),
    recordEvidence: protectedProcedure
      .input(z.object({ runId: z.string().min(1), evidenceClass: z.string().min(1), sourceHash: z.string().min(1), payload: z.unknown(), provenance: z.unknown() }))
      .mutation(({ ctx, input }) => recordKernelEvidence({ userId: ctx.user.id, ...input })),
    recordDrift: protectedProcedure
      .input(z.object({ runId: z.string().optional(), artifactSha256: z.string().min(1), baselineSha256: z.string().min(1), impact: z.array(z.string()) }))
      .mutation(({ ctx, input }) => recordKernelDrift({ userId: ctx.user.id, ...input })),
    registerEngine: protectedProcedure
      .input(z.object({ engineId: z.string().min(1), name: z.string().min(1), version: z.string().min(1), adapterKind: z.string().min(1), capabilities: z.array(z.string()) }))
      .mutation(({ ctx, input }) => registerKernelEngine({ userId: ctx.user.id, ...input })),
    registerTest: protectedProcedure
      .input(z.object({ testId: z.string().min(1), version: z.string().min(1), purpose: z.string().min(1), preconditions: z.array(z.string()), inputs: z.array(z.string()), inputHashes: z.array(z.string()), engineId: z.string().min(1), engineVersion: z.string().min(1), command: z.string().min(1), checks: z.array(z.string()), acceptance: z.array(z.string()), outputs: z.array(z.string()), evidenceClass: z.string().min(1), gatePolicy: z.string().min(1), protectedArtifact: z.string().optional() }))
      .mutation(({ ctx, input }) => registerKernelTest({ userId: ctx.user.id, ...input })),
  }),

  // TODO: add feature routers here, e.g.
  // todo: router({
  //   list: protectedProcedure.query(({ ctx }) =>
  //     db.getUserTodos(ctx.user.id)
  //   ),
  // }),
});

export type AppRouter = typeof appRouter;
