package squish

import (
	"fmt"
	"log/slog"
	"os"
	"runtime/pprof"
	"strings"
)

func startProfiling(cpuFilename string, memFilename string) (*os.File, *os.File, error) {
	cpuFile, err := os.Create(cpuFilename)
	if err != nil {
		return nil, nil, fmt.Errorf("create CPU profile %q: %w", cpuFilename, err)
	}
	slog.Debug("CPU profile will be saved", "path", cpuFilename)
	if err := pprof.StartCPUProfile(cpuFile); err != nil {
		cpuFile.Close()
		return nil, nil, fmt.Errorf("start CPU profile: %w", err)
	}

	memFile, err := os.Create(memFilename)
	if err != nil {
		pprof.StopCPUProfile()
		cpuFile.Close()
		return nil, nil, fmt.Errorf("create memory profile %q: %w", memFilename, err)
	}
	slog.Debug("Memory profile will be saved", "path", memFilename)
	return cpuFile, memFile, nil
}

func stopProfiling(cpuFile *os.File, memFile *os.File) error {
	pprof.StopCPUProfile()
	var errs []string
	if err := pprof.WriteHeapProfile(memFile); err != nil {
		errs = append(errs, fmt.Sprintf("write heap profile: %v", err))
	}
	if err := cpuFile.Close(); err != nil {
		errs = append(errs, fmt.Sprintf("close CPU profile: %v", err))
	}
	if err := memFile.Close(); err != nil {
		errs = append(errs, fmt.Sprintf("close memory profile: %v", err))
	}
	if len(errs) > 0 {
		return fmt.Errorf(strings.Join(errs, "; "))
	}
	return nil
}
